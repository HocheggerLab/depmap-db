"""Processor for phase-1 PRISM primary repurposing data."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd

from ...utils.constants import DEPMAP_FILES
from .base import BaseProcessor, ProcessingResult
from .prism_common import PrismCompoundMixin


class PrismPrimaryWideProcessor(PrismCompoundMixin, BaseProcessor):
    """Load the PRISM primary matrix as-published (compound x model).

    Phase 1 preserves the source orientation rather than rotating it into a long table.
    Compound metadata is loaded from the companion compound list when available.
    """

    DATASET_NAME = "PRISMPrimaryRepurposingExtended"
    TABLE_NAME = "drug_response_primary_wide"
    COMPOUND_LIST_DATASET = "PRISMPrimaryRepurposingCompoundList"

    def __init__(self) -> None:
        super().__init__(self.DATASET_NAME, batch_size=None)
        self.dataset_info = DEPMAP_FILES[self.DATASET_NAME]
        self.column_mapping: dict[str, str] = {}

    def get_table_name(self) -> str:
        return self.TABLE_NAME

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        warnings: list[str] = []
        if df.empty:
            raise ValueError("PRISM primary matrix is empty")
        if len(df.columns) < 2:
            raise ValueError(
                "PRISM primary matrix must contain one compound column plus model columns"
            )

        first_column = df.columns[0]
        missing_compounds = int(df[first_column].isna().sum())
        if missing_compounds > 0:
            warnings.append(
                f"Dropping {missing_compounds} rows with missing compound identifiers"
            )
            df = df.dropna(subset=[first_column]).copy()

        duplicate_compounds = int(df[first_column].duplicated().sum())
        if duplicate_compounds > 0:
            warnings.append(
                f"Found {duplicate_compounds} duplicate compound rows; keeping first occurrence"
            )
            df = df.drop_duplicates(subset=[first_column], keep="first").copy()

        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        return df

    def process_file(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        start_time = pd.Timestamp.now()

        try:
            df = pd.read_csv(file_path, low_memory=False)
            df, warnings = self.validate_data(df)
            compound_column = df.columns[0]
            model_columns = [str(column) for column in df.columns[1:]]

            # Clear wide table first to avoid FK conflicts when replacing the screen row.
            self.db_manager.execute(f"DROP TABLE IF EXISTS {self.TABLE_NAME}")

            self._ensure_screen_row(
                screen_id="PRISM_PRIMARY_24Q2",
                dataset_name=self.DATASET_NAME,
                source_filename=self.dataset_info.filename,
                screen_kind="primary",
                release_label=self.dataset_info.release_label_override
                or self.settings.depmap.release_label,
            )

            compound_metadata, metadata_warnings = (
                self._load_compound_metadata(file_path)
            )
            warnings.extend(metadata_warnings)
            self._replace_compounds_table(compound_metadata)
            target_df = self._resolve_compound_targets(
                compound_metadata,
                target_column="target_text",
                source_dataset=self.COMPOUND_LIST_DATASET,
                source_filename=DEPMAP_FILES[
                    self.COMPOUND_LIST_DATASET
                ].filename,
            )
            self._replace_compound_targets(
                target_df,
                source_dataset=self.COMPOUND_LIST_DATASET,
            )

            # Backfill minimal compound rows for any broad_ids in the matrix
            # that are missing from the companion compound list.
            matrix_broad_ids = set(df[compound_column].astype(str).unique())
            existing_rows = self.db_manager.execute(
                "SELECT broad_id FROM compounds"
            ).fetchall()
            existing_broad_ids = {row[0] for row in existing_rows}
            missing_ids = matrix_broad_ids - existing_broad_ids
            if missing_ids:
                warnings.append(
                    f"Backfilling {len(missing_ids)} compound rows not in compound list"
                )
                missing_df = pd.DataFrame(
                    {
                        "broad_id": list(missing_ids),
                        "source_dataset": self.DATASET_NAME,
                        "source_filename": self.dataset_info.filename,
                        "release_label": self.dataset_info.release_label_override,
                        "created_at": pd.Timestamp.now(),
                    }
                )
                for col in [
                    "compound_name",
                    "compound_synonyms",
                    "moa",
                    "target_text",
                    "smiles",
                    "phase",
                    "primary_screen_id",
                    "primary_dose_um",
                    "secondary_screen_id",
                    "secondary_row_name",
                    "secondary_passed_str_profiling",
                    "disease_area",
                    "indication",
                ]:
                    missing_df[col] = None
                compound_columns = [
                    "broad_id",
                    "compound_name",
                    "compound_synonyms",
                    "moa",
                    "target_text",
                    "smiles",
                    "phase",
                    "primary_screen_id",
                    "primary_dose_um",
                    "secondary_screen_id",
                    "secondary_row_name",
                    "secondary_passed_str_profiling",
                    "disease_area",
                    "indication",
                    "source_dataset",
                    "source_filename",
                    "release_label",
                    "created_at",
                ]
                missing_df = missing_df[compound_columns]
                temp_name = "prism_missing_compounds"
                conn = self.db_manager.connect()
                self.db_manager.execute(f"DROP VIEW IF EXISTS {temp_name}")
                conn.register(temp_name, missing_df)
                try:
                    cols_sql = ", ".join(compound_columns)
                    self.db_manager.execute(
                        f"INSERT INTO compounds ({cols_sql}) SELECT {cols_sql} FROM {temp_name}"
                    )
                finally:
                    conn.unregister(temp_name)

            self.create_wide_table_schema(model_columns)
            df_insert = df.rename(
                columns={compound_column: "broad_id", **self.column_mapping}
            ).copy()
            df_insert.insert(1, "screen_id", "PRISM_PRIMARY_24Q2")
            df_insert.insert(2, "created_at", pd.Timestamp.now())
            ordered_columns = [
                "broad_id",
                "screen_id",
                "created_at",
                *self.column_mapping.values(),
            ]
            df_insert = df_insert[ordered_columns]

            with NamedTemporaryFile(
                suffix=".parquet", delete=False
            ) as temp_file:
                temp_path = Path(temp_file.name)
            try:
                df_insert.to_parquet(temp_path, index=False)
                self.db_manager.execute(f"DELETE FROM {self.TABLE_NAME}")
                self.db_manager.execute(
                    f"""
                    INSERT INTO {self.TABLE_NAME}
                    SELECT * FROM read_parquet('{temp_path}')
                    """
                )
            finally:
                temp_path.unlink(missing_ok=True)

            self._record_import(file_path, len(df_insert), "success")
            processing_time = (pd.Timestamp.now() - start_time).total_seconds()
            return ProcessingResult(
                processor_name=self.__class__.__name__,
                dataset_name=self.dataset_name,
                source_file=file_path,
                records_processed=len(df_insert),
                records_inserted=len(df_insert),
                records_updated=0,
                records_skipped=0,
                processing_time_seconds=processing_time,
                status="success",
                warnings=warnings,
            )
        except (OSError, ValueError, RuntimeError) as e:
            processing_time = (pd.Timestamp.now() - start_time).total_seconds()
            error_msg = f"Processing failed: {e}"
            self.logger.error("%s", error_msg)
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

    def create_wide_table_schema(self, model_ids: list[str]) -> None:
        self.db_manager.execute(f"DROP TABLE IF EXISTS {self.TABLE_NAME}")
        self.column_mapping = {
            model_id: self._sanitize_model_column(model_id)
            for model_id in model_ids
        }
        model_columns_sql = ",\n                ".join(
            f'"{column_name}" DOUBLE'
            for column_name in self.column_mapping.values()
        )
        self.db_manager.execute(
            f"""
            CREATE TABLE {self.TABLE_NAME} (
                broad_id VARCHAR PRIMARY KEY,
                screen_id VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                {model_columns_sql},
                FOREIGN KEY (broad_id) REFERENCES compounds(broad_id),
                FOREIGN KEY (screen_id) REFERENCES drug_screens(screen_id)
            )
            """
        )

    def _sanitize_model_column(self, model_id: str) -> str:
        return f"model_{str(model_id).replace('-', '_').replace('.', '_')}"

    def _load_compound_metadata(
        self, matrix_path: Path
    ) -> tuple[pd.DataFrame, list[str]]:
        warnings: list[str] = []
        compound_list_info = DEPMAP_FILES[self.COMPOUND_LIST_DATASET]
        compound_list_path = matrix_path.with_name(compound_list_info.filename)

        if not compound_list_path.exists():
            warnings.append(
                "Primary compound list not found alongside matrix; creating minimal compounds rows from broad IDs only"
            )
            matrix_df = pd.read_csv(matrix_path, usecols=[0], low_memory=False)
            first_column = matrix_df.columns[0]
            compounds_df = pd.DataFrame(
                {
                    "broad_id": matrix_df[first_column].astype(str),
                    "compound_name": None,
                    "compound_synonyms": None,
                    "moa": None,
                    "target_text": None,
                    "smiles": None,
                    "phase": None,
                    "primary_screen_id": None,
                    "primary_dose_um": None,
                    "secondary_screen_id": None,
                    "secondary_row_name": None,
                    "secondary_passed_str_profiling": None,
                    "disease_area": None,
                    "indication": None,
                    "source_dataset": self.DATASET_NAME,
                    "source_filename": self.dataset_info.filename,
                    "release_label": self.dataset_info.release_label_override,
                }
            )
            return compounds_df.drop_duplicates(subset=["broad_id"]), warnings

        raw = pd.read_csv(compound_list_path, low_memory=False)
        required_columns = {"IDs", "Drug.Name"}
        missing = required_columns - set(raw.columns)
        if missing:
            raise ValueError(
                f"Primary compound list missing required columns: {', '.join(sorted(missing))}"
            )

        compounds_df = raw.rename(
            columns={
                "IDs": "broad_id",
                "Drug.Name": "compound_name",
                "Synonyms": "compound_synonyms",
                "MOA": "moa",
                "repurposing_target": "target_text",
                "screen": "primary_screen_id",
                "dose": "primary_dose_um",
            }
        )[
            [
                "broad_id",
                "compound_name",
                "compound_synonyms",
                "moa",
                "target_text",
                "primary_screen_id",
                "primary_dose_um",
            ]
        ].copy()
        compounds_df["smiles"] = None
        compounds_df["phase"] = None
        compounds_df["secondary_screen_id"] = None
        compounds_df["secondary_row_name"] = None
        compounds_df["secondary_passed_str_profiling"] = None
        compounds_df["disease_area"] = None
        compounds_df["indication"] = None
        compounds_df["source_dataset"] = self.COMPOUND_LIST_DATASET
        compounds_df["source_filename"] = compound_list_info.filename
        compounds_df["release_label"] = (
            compound_list_info.release_label_override
        )
        compounds_df = compounds_df.dropna(
            subset=["broad_id"]
        ).drop_duplicates(subset=["broad_id"], keep="first")
        return compounds_df, warnings

    def _ensure_screen_row(
        self,
        *,
        screen_id: str,
        dataset_name: str,
        source_filename: str,
        screen_kind: str,
        release_label: str | None,
    ) -> None:
        self.db_manager.execute(
            "DELETE FROM drug_screens WHERE screen_id = ?",
            [screen_id],
        )
        self.db_manager.execute(
            """
            INSERT INTO drug_screens (
                screen_id,
                screen_kind,
                dataset_name,
                source_filename,
                release_label,
                release_track,
                default_secondary_summary_metric,
                notes
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                screen_id,
                screen_kind,
                dataset_name,
                source_filename,
                release_label,
                self.dataset_info.release_track,
                "auc",
                (
                    "Primary PRISM matrix stored as published in wide compound-by-model orientation. "
                    "Secondary analyses should generally prefer AUC when using fitted dose-response summaries."
                ),
            ],
        )
