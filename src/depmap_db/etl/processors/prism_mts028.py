"""Processor for the internal MTS028 HELFRID HOCHEGGER PRISM-like c-604 screen."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd

from ...utils.constants import DatasetInfo
from .base import BaseProcessor, ProcessingResult
from .prism_secondary import PrismSecondaryProcessor


class PrismMTS028Processor(PrismSecondaryProcessor):
    """Load the internal MTS028 c-604 screen into the canonical PRISM tables.

    Design choice:
    - reuse ``drug_response_secondary`` for fitted summary rows
    - add ``drug_response_secondary_dose`` for collapsed dose-level l2fc values
    - reuse the shared compounds/drug_screens abstractions
    """

    DATASET_NAME = "PRISMCustomMTS028HELFRIDHOCHEGGER"
    DOSE_TABLE_NAME = "drug_response_secondary_dose"
    COMPOUND_MOA = "Greatwall/MASTL kinase inhibitor"
    COMPOUND_TARGET_TEXT = "MASTL"

    def __init__(self) -> None:
        BaseProcessor.__init__(self, self.DATASET_NAME, batch_size=None)
        self.dataset_info = DatasetInfo(
            filename="MTS028_HELFRID_HOCHEGGER_DRC_TABLE.csv",
            description=(
                "Internal PRISM-like dose-response screen for c-604 with DRC "
                "summaries and collapsed dose-level l2fc values."
            ),
            table_name=self.TABLE_NAME,
            priority=2,
            modality="drug_sensitivity",
            release_track="prism_custom_internal",
            release_label_override="Internal PRISM-like c-604 MTS028 screen",
        )

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        warnings: list[str] = []
        required = {
            "pert_id",
            "pert_name",
            "depmap_id",
            "screen",
            "auc",
            "ic50",
            "successful_fit",
            "x_project_id",
            "minimum_dose",
            "maximum_dose",
        }
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                "MTS028 DRC table missing required columns: "
                + ", ".join(sorted(missing))
            )
        return df, warnings

    def process_file(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        start_time = pd.Timestamp.now()

        try:
            df = pd.read_csv(file_path, low_memory=False)
            df, warnings = self.validate_data(df)

            self._ensure_custom_schema_compatibility()
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

            summary_df = self._build_secondary_response_dataframe(df)
            dose_df, dose_warnings = self._load_collapsed_dose_dataframe(file_path)
            warnings.extend(dose_warnings)

            response_model_ids = set(
                summary_df["model_id"].dropna().astype(str).unique()
            )
            existing_model_rows = self.db_manager.execute(
                "SELECT model_id FROM models"
            ).fetchall()
            existing_model_ids = {row[0] for row in existing_model_rows}
            missing_model_ids = response_model_ids - existing_model_ids
            if missing_model_ids:
                warnings.append(
                    f"Backfilling {len(missing_model_ids)} stub model rows from MTS028 custom screen"
                )
                for mid in missing_model_ids:
                    self.db_manager.execute(
                        "INSERT INTO models (model_id, cell_line_name) VALUES (?, ?)",
                        [mid, mid],
                    )

            response_broad_ids = set(
                summary_df["broad_id"].dropna().astype(str).unique()
            )
            existing_compound_rows = self.db_manager.execute(
                "SELECT broad_id FROM compounds"
            ).fetchall()
            existing_compound_ids = {row[0] for row in existing_compound_rows}
            missing_broad_ids = response_broad_ids - existing_compound_ids
            if missing_broad_ids:
                warnings.append(
                    f"Backfilling {len(missing_broad_ids)} stub compound rows from MTS028 custom screen"
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

            self._replace_secondary_table(summary_df)
            self._replace_secondary_dose_table(dose_df)

            self._record_import(file_path, len(summary_df), "success")
            processing_time = (pd.Timestamp.now() - start_time).total_seconds()
            return ProcessingResult(
                processor_name=self.__class__.__name__,
                dataset_name=self.dataset_name,
                source_file=file_path,
                records_processed=len(summary_df),
                records_inserted=len(summary_df),
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

    def _ensure_custom_schema_compatibility(self) -> None:
        secondary_alters = [
            "ALTER TABLE drug_response_secondary ADD COLUMN IF NOT EXISTS fit_name VARCHAR",
            "ALTER TABLE drug_response_secondary ADD COLUMN IF NOT EXISTS successful_fit BOOLEAN",
            "ALTER TABLE drug_response_secondary ADD COLUMN IF NOT EXISTS auc_riemann DOUBLE",
            "ALTER TABLE drug_response_secondary ADD COLUMN IF NOT EXISTS minimum_dose_um DOUBLE",
            "ALTER TABLE drug_response_secondary ADD COLUMN IF NOT EXISTS maximum_dose_um DOUBLE",
            "ALTER TABLE drug_response_secondary ADD COLUMN IF NOT EXISTS source_project_id VARCHAR",
        ]
        for sql in secondary_alters:
            self.db_manager.execute(sql)

        self.db_manager.execute(
            """
            CREATE TABLE IF NOT EXISTS drug_response_secondary_dose (
                dose_response_id VARCHAR PRIMARY KEY,
                response_id VARCHAR NOT NULL,
                broad_id VARCHAR NOT NULL,
                model_id VARCHAR NOT NULL,
                screen_id VARCHAR NOT NULL,
                dose_um DOUBLE NOT NULL,
                dose_unit VARCHAR,
                median_l2fc DOUBLE,
                median_l2fc_uncorrected DOUBLE,
                num_bio_reps INTEGER,
                pool_id VARCHAR,
                lua VARCHAR,
                cell_set VARCHAR,
                growth_pattern VARCHAR,
                pert_type VARCHAR,
                pert_vehicle VARCHAR,
                pert_plate VARCHAR,
                day INTEGER,
                source_project_id VARCHAR,
                source_dataset VARCHAR NOT NULL,
                source_filename VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (response_id) REFERENCES drug_response_secondary(response_id),
                FOREIGN KEY (broad_id) REFERENCES compounds(broad_id),
                FOREIGN KEY (model_id) REFERENCES models(model_id),
                FOREIGN KEY (screen_id) REFERENCES drug_screens(screen_id)
            )
            """
        )

    def _build_compounds_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        compounds_df = pd.DataFrame(
            {
                "broad_id": df["pert_id"],
                "compound_name": df["pert_name"],
                "compound_synonyms": None,
                "moa": self.COMPOUND_MOA,
                "target_text": self.COMPOUND_TARGET_TEXT,
                "smiles": None,
                "phase": None,
                "primary_screen_id": None,
                "primary_dose_um": None,
                "secondary_screen_id": df["screen"],
                "secondary_row_name": None,
                "secondary_passed_str_profiling": None,
                "disease_area": None,
                "indication": None,
                "source_dataset": self.DATASET_NAME,
                "source_filename": self.dataset_info.filename,
                "release_label": self.dataset_info.release_label_override,
            }
        )
        return compounds_df.dropna(subset=["broad_id"]).drop_duplicates(
            subset=["broad_id"], keep="first"
        )

    def _build_secondary_response_dataframe(
        self, df: pd.DataFrame
    ) -> pd.DataFrame:
        clean = df.dropna(subset=["depmap_id", "pert_id", "screen"]).copy()
        clean["response_id"] = (
            clean["pert_id"].astype(str)
            + "::"
            + clean["depmap_id"].astype(str)
            + "::"
            + clean["screen"].astype(str)
        )
        clean = clean.drop_duplicates(subset=["response_id"], keep="first")
        clean["created_at"] = pd.Timestamp.now()
        clean["source_dataset"] = self.DATASET_NAME
        clean["source_filename"] = self.dataset_info.filename
        clean["ec50"] = None
        clean["passed_str_profiling"] = None
        clean["row_name"] = None
        clean["compound_name"] = clean["pert_name"]
        clean["moa"] = self.COMPOUND_MOA
        clean["target_text"] = self.COMPOUND_TARGET_TEXT
        clean["disease_area"] = None
        clean["indication"] = None
        clean["smiles"] = None
        clean["phase"] = None

        ordered = pd.DataFrame(
            {
                "response_id": clean["response_id"],
                "broad_id": clean["pert_id"],
                "model_id": clean["depmap_id"],
                "ccle_name": None,
                "screen_id": clean["screen"],
                "upper_limit": clean["upper_limit"],
                "lower_limit": clean["lower_limit"],
                "slope": clean["slope"],
                "r2": clean["frac_var_explained"],
                "auc": clean["auc"],
                "ec50": clean["ec50"],
                "ic50": clean["ic50"],
                "fit_name": clean["fit_name"],
                "successful_fit": clean["successful_fit"],
                "auc_riemann": clean["auc_riemann"],
                "minimum_dose_um": clean["minimum_dose"],
                "maximum_dose_um": clean["maximum_dose"],
                "source_project_id": clean["x_project_id"],
                "passed_str_profiling": clean["passed_str_profiling"],
                "row_name": clean["row_name"],
                "compound_name": clean["compound_name"],
                "moa": clean["moa"],
                "target_text": clean["target_text"],
                "disease_area": clean["disease_area"],
                "indication": clean["indication"],
                "smiles": clean["smiles"],
                "phase": clean["phase"],
                "source_dataset": clean["source_dataset"],
                "source_filename": clean["source_filename"],
                "created_at": clean["created_at"],
            }
        )
        return ordered

    def _ensure_screen_rows(self, df: pd.DataFrame) -> None:
        release_label = self.dataset_info.release_label_override
        screens = sorted(df["screen"].dropna().astype(str).unique().tolist())
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
                    "secondary_custom",
                    self.DATASET_NAME,
                    self.dataset_info.filename,
                    release_label,
                    self.dataset_info.release_track,
                    "auc",
                    (
                        "Internal PRISM-like c-604 dose-response screen. "
                        "Fitted summaries are stored in drug_response_secondary and "
                        "collapsed dose-level l2fc values are stored in "
                        "drug_response_secondary_dose."
                    ),
                ],
            )

    def _load_collapsed_dose_dataframe(
        self, drc_path: Path
    ) -> tuple[pd.DataFrame, list[str]]:
        warnings: list[str] = []
        collapsed_path = (
            drc_path.parent.parent
            / "build"
            / "MTS028_HELFRID_HOCHEGGER_collapsed_l2fc.csv"
        )
        if not collapsed_path.exists():
            warnings.append(
                "Collapsed dose-level file not found; dose-level table was not populated"
            )
            return pd.DataFrame(), warnings

        dose_df = pd.read_csv(collapsed_path, low_memory=False)
        required = {
            "pert_id",
            "depmap_id",
            "screen",
            "pert_dose",
            "pert_dose_unit",
            "median_l2fc",
            "median_l2fc_uncorrected",
            "num_bio_reps",
            "x_project_id",
        }
        missing = required - set(dose_df.columns)
        if missing:
            raise ValueError(
                "Collapsed dose file missing required columns: "
                + ", ".join(sorted(missing))
            )

        clean = dose_df.dropna(
            subset=["pert_id", "depmap_id", "screen", "pert_dose"]
        ).copy()
        clean["response_id"] = (
            clean["pert_id"].astype(str)
            + "::"
            + clean["depmap_id"].astype(str)
            + "::"
            + clean["screen"].astype(str)
        )
        clean["dose_response_id"] = (
            clean["response_id"]
            + "::"
            + clean["pert_dose"].astype(str)
        )
        clean = clean.drop_duplicates(
            subset=["dose_response_id"], keep="first"
        )
        clean["source_dataset"] = self.DATASET_NAME
        clean["source_filename"] = collapsed_path.name
        clean["created_at"] = pd.Timestamp.now()

        ordered = pd.DataFrame(
            {
                "dose_response_id": clean["dose_response_id"],
                "response_id": clean["response_id"],
                "broad_id": clean["pert_id"],
                "model_id": clean["depmap_id"],
                "screen_id": clean["screen"],
                "dose_um": clean["pert_dose"],
                "dose_unit": clean["pert_dose_unit"],
                "median_l2fc": clean["median_l2fc"],
                "median_l2fc_uncorrected": clean[
                    "median_l2fc_uncorrected"
                ],
                "num_bio_reps": clean["num_bio_reps"],
                "pool_id": clean.get("pool_id"),
                "lua": clean.get("lua"),
                "cell_set": clean.get("cell_set"),
                "growth_pattern": clean.get("growth_pattern"),
                "pert_type": clean.get("pert_type"),
                "pert_vehicle": clean.get("pert_vehicle"),
                "pert_plate": clean.get("pert_plate"),
                "day": clean.get("day"),
                "source_project_id": clean["x_project_id"],
                "source_dataset": clean["source_dataset"],
                "source_filename": clean["source_filename"],
                "created_at": clean["created_at"],
            }
        )
        return ordered, warnings

    def _replace_secondary_table(self, df_insert: pd.DataFrame) -> None:
        columns = ", ".join(df_insert.columns)
        with NamedTemporaryFile(suffix=".parquet", delete=False) as temp_file:
            temp_path = Path(temp_file.name)
        try:
            df_insert.to_parquet(temp_path, index=False)
            self.db_manager.execute(
                "DELETE FROM drug_response_secondary WHERE source_dataset = ?",
                [self.DATASET_NAME],
            )
            self.db_manager.execute(
                f"""
                INSERT INTO {self.TABLE_NAME} ({columns})
                SELECT {columns} FROM read_parquet('{temp_path}')
                """
            )
        finally:
            temp_path.unlink(missing_ok=True)

    def _replace_secondary_dose_table(self, dose_df: pd.DataFrame) -> None:
        self.db_manager.execute(
            "DELETE FROM drug_response_secondary_dose WHERE source_dataset = ?",
            [self.DATASET_NAME],
        )
        if dose_df.empty:
            return

        columns = ", ".join(dose_df.columns)
        with NamedTemporaryFile(suffix=".parquet", delete=False) as temp_file:
            temp_path = Path(temp_file.name)
        try:
            dose_df.to_parquet(temp_path, index=False)
            self.db_manager.execute(
                f"""
                INSERT INTO {self.DOSE_TABLE_NAME} ({columns})
                SELECT {columns} FROM read_parquet('{temp_path}')
                """
            )
        finally:
            temp_path.unlink(missing_ok=True)
