"""Base processor class for ETL operations."""

import contextlib
import uuid
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pandas as pd

from ...config import get_logger, get_settings
from ...database import get_db_manager

logger = get_logger(__name__)


@dataclass
class ProcessingResult:
    """Result of a data processing operation."""

    processor_name: str
    dataset_name: str
    source_file: Path
    records_processed: int
    records_inserted: int
    records_updated: int
    records_skipped: int
    processing_time_seconds: float
    status: str  # 'success', 'partial', 'failed'
    error_message: str | None = None
    warnings: list[str] | None = None

    def __post_init__(self) -> None:
        if self.warnings is None:
            self.warnings = []


class BaseProcessor(ABC):
    """Abstract base class for data processors."""

    def __init__(self, dataset_name: str, batch_size: int | None = None):
        """Initialize processor.

        Args:
            dataset_name: Name of the dataset being processed
            batch_size: Optional batch size for processing. Uses settings default if None.
        """
        self.dataset_name = dataset_name
        self.settings = get_settings()
        self.batch_size = batch_size or self.settings.depmap.batch_size
        self.db_manager = get_db_manager()
        self.logger = get_logger(self.__class__.__name__)

    @abstractmethod
    def get_table_name(self) -> str:
        """Get the target database table name."""

    @abstractmethod
    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Validate and clean the input data.

        Args:
            df: Input DataFrame

        Returns:
            Tuple of (cleaned_dataframe, warnings_list)
        """

    @abstractmethod
    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform data to match database schema.

        Args:
            df: Input DataFrame

        Returns:
            Transformed DataFrame ready for database insertion
        """

    def load_file(self, file_path: Path) -> pd.DataFrame:
        """Load data from file.

        Args:
            file_path: Path to the data file

        Returns:
            DataFrame with loaded data
        """
        self.logger.info("Loading data from %s", file_path)

        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")

        # Detect file type and load accordingly
        if file_path.suffix.lower() == ".csv":
            df = pd.read_csv(file_path, low_memory=False)
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")

        self.logger.info("Loaded %s records from %s", len(df), file_path.name)
        return df

    def get_existing_records(self) -> set[str] | set[tuple[str, str]]:
        """Get set of existing record identifiers to avoid duplicates.

        Returns:
            Set of existing record identifiers
        """
        self.get_table_name()

        try:
            # This will be overridden by subclasses with appropriate logic
            return set()

        except RuntimeError as e:
            self.logger.warning("Could not retrieve existing records: %s", e)
            return set()

    def insert_batch(self, df: pd.DataFrame) -> tuple[int, int]:
        """Insert a batch of records into the database.

        Args:
            df: DataFrame to insert

        Returns:
            Tuple of (records_inserted, records_updated)
        """
        if df.empty:
            return 0, 0

        table_name = self.get_table_name()

        try:
            # Use DuckDB's efficient CSV import when possible
            # For now, we'll use pandas to_sql equivalent

            # Convert DataFrame to records for insertion
            records = df.to_dict("records")

            if not records:
                return 0, 0

            # Build parameterized insert query
            columns = list(df.columns)
            placeholders = ", ".join(["?" for _ in columns])
            column_names = ", ".join(columns)

            insert_query = f"""
            INSERT OR REPLACE INTO {table_name} ({column_names})
            VALUES ({placeholders})
            """

            # Convert records to parameter tuples
            parameters = [
                tuple(record[col] for col in columns) for record in records
            ]

            # Execute batch insert
            self.db_manager.execute_many(insert_query, parameters)

            return len(records), 0  # For now, treating all as inserts

        except RuntimeError as e:
            self.logger.error("Failed to insert batch: %s", e)
            raise

    def process_file(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        """Process a data file from start to finish.

        Args:
            file_path: Path to the data file
            force_reload: If True, reload data even if already processed

        Returns:
            ProcessingResult with operation details
        """
        start_time = datetime.now()

        try:
            # Check if already processed (unless force_reload)
            if not force_reload:
                existing_count = self._get_existing_record_count()
                if existing_count > 0:
                    self.logger.info(
                        "Table %s already has %s records",
                        self.get_table_name(),
                        existing_count,
                    )
                    if not self._should_reprocess():
                        return ProcessingResult(
                            processor_name=self.__class__.__name__,
                            dataset_name=self.dataset_name,
                            source_file=file_path,
                            records_processed=0,
                            records_inserted=0,
                            records_updated=0,
                            records_skipped=existing_count,
                            processing_time_seconds=0.0,
                            status="skipped",
                        )

            # Load data
            df = self.load_file(file_path)
            total_records = len(df)

            # Validate data
            df, warnings = self.validate_data(df)

            # Transform data
            df = self.transform_data(df)

            # Process in batches
            records_inserted = 0
            records_updated = 0

            for i in range(0, len(df), self.batch_size):
                batch_df = df.iloc[i : i + self.batch_size]
                batch_inserted, batch_updated = self.insert_batch(batch_df)
                records_inserted += batch_inserted
                records_updated += batch_updated

                if i + self.batch_size < len(df):
                    self.logger.info(
                        "Processed %s/%s records",
                        i + len(batch_df),
                        len(df),
                    )

            # Record the import in data_imports table
            self._record_import(
                file_path, records_inserted + records_updated, "success"
            )

            processing_time = (datetime.now() - start_time).total_seconds()

            return ProcessingResult(
                processor_name=self.__class__.__name__,
                dataset_name=self.dataset_name,
                source_file=file_path,
                records_processed=total_records,
                records_inserted=records_inserted,
                records_updated=records_updated,
                records_skipped=total_records - len(df),
                processing_time_seconds=processing_time,
                status="success",
                warnings=warnings,
            )

        except (ValueError, OSError, RuntimeError) as e:
            processing_time = (datetime.now() - start_time).total_seconds()
            error_msg = f"Processing failed: {e}"
            self.logger.error("%s", error_msg)

            # Record failed import
            with contextlib.suppress(RuntimeError):
                self._record_import(file_path, 0, "failed", error_msg)

            return ProcessingResult(
                processor_name=self.__class__.__name__,
                dataset_name=self.dataset_name,
                source_file=file_path,
                records_processed=0,
                records_inserted=0,
                records_updated=0,
                records_skipped=0,
                processing_time_seconds=processing_time,
                status="failed",
                error_message=error_msg,
            )

    def _get_existing_record_count(self) -> int:
        """Get count of existing records in target table."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return 0

            result = self.db_manager.execute(
                f"SELECT COUNT(*) FROM {table_name}"
            )
            row = result.fetchone()
            return int(row[0]) if row else 0

        except RuntimeError:
            return 0

    def _should_reprocess(self) -> bool:
        """Determine if data should be reprocessed even when records exist.

        Subclasses can override this for custom logic.
        """
        return False

    def _record_import(
        self,
        file_path: Path,
        records_imported: int,
        status: str,
        error_message: str | None = None,
    ) -> None:
        """Record import operation in data_imports table."""
        import_id = str(uuid.uuid4())

        try:
            self.db_manager.execute(
                """
                INSERT INTO data_imports (
                    import_id, import_type, source_file, records_imported,
                    status, started_at, completed_at, error_message
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
                [
                    import_id,
                    self.dataset_name,
                    str(file_path),
                    records_imported,
                    status,
                    datetime.now().isoformat(),
                    datetime.now().isoformat(),
                    error_message,
                ],
            )

        except RuntimeError as e:
            self.logger.error("Failed to record import: %s", e)

    def clear_table(self) -> int:
        """Clear all data from the target table.

        Returns:
            Number of records deleted
        """
        table_name = self.get_table_name()

        try:
            # Get count before deletion
            count_result = self.db_manager.execute(
                f"SELECT COUNT(*) FROM {table_name}"
            )
            row = count_result.fetchone()
            record_count: int = int(row[0]) if row else 0

            # Delete all records
            self.db_manager.execute(f"DELETE FROM {table_name}")

            self.logger.info(
                "Cleared %s records from %s", record_count, table_name
            )
            return record_count

        except RuntimeError as e:
            self.logger.error("Failed to clear table %s: %s", table_name, e)
            raise
