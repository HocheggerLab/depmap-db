"""Processor for phase-1 PRISM secondary dose-response data."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd

from ...utils.constants import DEPMAP_FILES
from .base import BaseProcessor, ProcessingResult
from .prism_common import PrismCompoundMixin


class PrismSecondaryProcessor(PrismCompoundMixin, BaseProcessor):
    """Load the PRISM secondary dose-response parameter table.

    Phase 1 stores the long table canonically, preserving AUC/EC50/IC50 and fit
    terms. AUC is documented as the default summary for most first-pass analyses,
    but the full fitted-curve parameters remain available.
    """

    DATASET_NAME = "PRISMSecondaryDoseResponseCurveParameters"
    TABLE_NAME = "drug_response_secondary"

    def __init__(self) -> None:
        super().__init__(self.DATASET_NAME, batch_size=None)
        self.dataset_info = DEPMAP_FILES[self.DATASET_NAME]

    def get_table_name(self) -> str:
        return self.TABLE_NAME

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        warnings: list[str] = []
        required = {
            "broad_id",
            "depmap_id",
            "screen_id",
            "auc",
            "ec50",
            "ic50",
            "name",
            "target",
        }
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"PRISM secondary table missing required columns: {', '.join(sorted(missing))}"
            )
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

            self._ensure_screen_rows(df)
            compounds_df = self._build_compounds_dataframe(df)
            self._replace_compounds_table(compounds_df)
            target_df = self._resolve_compound_targets(
                compounds_df,
                target_column="target_text",
                source_dataset=self.DATASET_NAME,
                source_filename=self.dataset_info.filename,
            )
            self._replace_compound_targets(
                target_df, source_dataset=self.DATASET_NAME
            )

            df_insert = self._build_secondary_response_dataframe(df)

            # Backfill stub model rows for model_ids not in current models table.
            response_model_ids = set(
                df_insert["model_id"].dropna().astype(str).unique()
            )
            existing_model_rows = self.db_manager.execute(
                "SELECT model_id FROM models"
            ).fetchall()
            existing_model_ids = {row[0] for row in existing_model_rows}
            missing_model_ids = response_model_ids - existing_model_ids
            if missing_model_ids:
                warnings.append(
                    f"Backfilling {len(missing_model_ids)} stub model rows from PRISM secondary data"
                )
                for mid in missing_model_ids:
                    self.db_manager.execute(
                        "INSERT INTO models (model_id, cell_line_name) VALUES (?, ?)",
                        [mid, mid],
                    )

            # Backfill stub compound rows for broad_ids not in compounds table.
            response_broad_ids = set(
                df_insert["broad_id"].dropna().astype(str).unique()
            )
            existing_compound_rows = self.db_manager.execute(
                "SELECT broad_id FROM compounds"
            ).fetchall()
            existing_compound_ids = {row[0] for row in existing_compound_rows}
            missing_broad_ids = response_broad_ids - existing_compound_ids
            if missing_broad_ids:
                warnings.append(
                    f"Backfilling {len(missing_broad_ids)} stub compound rows from PRISM secondary data"
                )
                for bid in missing_broad_ids:
                    self.db_manager.execute(
                        "INSERT INTO compounds (broad_id, source_dataset, source_filename, release_label) VALUES (?, ?, ?, ?)",
                        [
                            bid,
                            self.DATASET_NAME,
                            self.dataset_info.filename,
                            self.dataset_info.release_label_override,
                        ],
                    )

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

    def _build_compounds_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        compounds_df = df.rename(
            columns={
                "name": "compound_name",
                "moa": "moa",
                "target": "target_text",
                "screen_id": "secondary_screen_id",
                "row_name": "secondary_row_name",
                "passed_str_profiling": "secondary_passed_str_profiling",
                "disease.area": "disease_area",
                "indication": "indication",
                "phase": "phase",
                "smiles": "smiles",
            }
        )[
            [
                "broad_id",
                "compound_name",
                "moa",
                "target_text",
                "smiles",
                "phase",
                "secondary_screen_id",
                "secondary_row_name",
                "secondary_passed_str_profiling",
                "disease_area",
                "indication",
            ]
        ].copy()
        compounds_df["compound_synonyms"] = None
        compounds_df["primary_screen_id"] = None
        compounds_df["primary_dose_um"] = None
        compounds_df["source_dataset"] = self.DATASET_NAME
        compounds_df["source_filename"] = self.dataset_info.filename
        compounds_df["release_label"] = (
            self.dataset_info.release_label_override
        )
        compounds_df = compounds_df.dropna(
            subset=["broad_id"]
        ).drop_duplicates(subset=["broad_id"], keep="first")
        return compounds_df[
            [
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
            ]
        ]

    def _build_secondary_response_dataframe(
        self, df: pd.DataFrame
    ) -> pd.DataFrame:
        clean = df.dropna(subset=["depmap_id", "broad_id"]).copy()
        clean["response_id"] = (
            clean["broad_id"].astype(str)
            + "::"
            + clean["depmap_id"].astype(str)
            + "::"
            + clean["screen_id"].astype(str)
        )
        clean = clean.drop_duplicates(subset=["response_id"], keep="first")
        clean = clean.rename(
            columns={
                "depmap_id": "model_id",
                "name": "compound_name",
                "moa": "moa",
                "target": "target_text",
                "disease.area": "disease_area",
                "passed_str_profiling": "passed_str_profiling",
            }
        )
        clean["created_at"] = pd.Timestamp.now()
        clean["source_dataset"] = self.DATASET_NAME
        clean["source_filename"] = self.dataset_info.filename
        clean["fit_name"] = None
        clean["successful_fit"] = None
        clean["auc_riemann"] = None
        clean["minimum_dose_um"] = None
        clean["maximum_dose_um"] = None
        clean["source_project_id"] = None
        ordered = [
            "response_id",
            "broad_id",
            "model_id",
            "ccle_name",
            "screen_id",
            "upper_limit",
            "lower_limit",
            "slope",
            "r2",
            "auc",
            "ec50",
            "ic50",
            "fit_name",
            "successful_fit",
            "auc_riemann",
            "minimum_dose_um",
            "maximum_dose_um",
            "source_project_id",
            "passed_str_profiling",
            "row_name",
            "compound_name",
            "moa",
            "target_text",
            "disease_area",
            "indication",
            "smiles",
            "phase",
            "source_dataset",
            "source_filename",
            "created_at",
        ]
        return clean[ordered]

    def _ensure_screen_rows(self, df: pd.DataFrame) -> None:
        release_label = (
            self.dataset_info.release_label_override
            or self.settings.depmap.release_label
        )
        screens = sorted(
            df["screen_id"].dropna().astype(str).unique().tolist()
        )
        for screen_id in screens:
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
                    "secondary",
                    self.DATASET_NAME,
                    self.dataset_info.filename,
                    release_label,
                    self.dataset_info.release_track,
                    "auc",
                    (
                        "Secondary PRISM dose-response fits retain auc/ec50/ic50 plus fit terms. "
                        "AUC is the recommended default summary for phase-1 analyses because it is "
                        "more stable than relying on a single fitted concentration summary alone."
                    ),
                ],
            )
