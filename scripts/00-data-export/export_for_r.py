# %%
"""Export the common per-cell / per-patient tables that every R script needs.

R has no ai4bmr-datasets binding, so every R script in this repo reads Parquet
files from `EXPORT_DIR` instead. This script produces the tables shared across
figure branches: `metadata.parquet` (per-cell labels), `clinical.parquet`
(per-ROI/patient clinical annotations), and `intensity_normalized.parquet`
(arcsinh + min-max normalized marker intensities, same transform used for
clustering). Figure-specific score tables (e.g. `scores.parquet`,
`survival-*.parquet`, niche abundance tables) are produced by each figure
branch's own scripts, not here.

Requires the full labeling pipeline to have already run (see
REPRODUCIBILITY.md) so that `01_raw/annotations/labels.parquet` exists.
"""
import os
from pathlib import Path

from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import prepare_data, resolve_base_dir

# Reported 3x consistently in the paper (Abstract, Results, Methods). If your
# materialized dataset diverges, that's worth investigating before trusting
# downstream figures — see the "Known discrepancies" section of
# REPRODUCIBILITY.md for a related count found in ai4bmr-datasets' own code.
EXPECTED_CELL_COUNT = 2_191_967
EXPECTED_PATIENT_COUNT = 195  # initial TMA cohort
EXPECTED_TUMOR_PATIENT_COUNT = 190  # final analytical cohort (is_tumor == "yes")


def main(base_dir: Path | None = None, export_dir: Path | None = None):
    from ai4bmr_datasets import PCa

    base_dir = Path(base_dir).expanduser() if base_dir else resolve_base_dir()
    export_dir = Path(export_dir).expanduser() if export_dir else Path(
        os.environ.get("EXPORT_DIR", base_dir / "0-export")
    )
    export_dir.mkdir(parents=True, exist_ok=True)

    # %% per-cell labels + per-ROI clinical annotations
    ds = PCa(
        base_dir=base_dir,
        image_version="filtered",
        mask_version="annotated",
        load_metadata=True,
        load_intensity=False,
        load_spatial=False,
    )
    ds.setup()

    metadata = ds.metadata.copy()
    clinical = ds.clinical.copy()

    assert len(metadata) == EXPECTED_CELL_COUNT, (
        f"cell count {len(metadata)} != paper-reported {EXPECTED_CELL_COUNT}; "
        "see REPRODUCIBILITY.md Known discrepancies before trusting downstream figures"
    )
    assert "pat_id" in clinical.columns, "clinical table is missing pat_id"

    # `clinical` (from ds.clinical) can contain ROIs with no processed/labeled
    # cells at all (e.g. failed segmentation) — restrict to ROIs actually
    # present in `metadata` before counting patients, which is what the paper's
    # 195 / 190 figures describe. See REPRODUCIBILITY.md Known discrepancies.
    roi_ids_with_cells = metadata.index.get_level_values("sample_id").unique()
    clinical_with_cells = clinical.loc[clinical.index.isin(roi_ids_with_cells)]

    n_patients = clinical_with_cells["pat_id"].nunique()
    assert n_patients == EXPECTED_PATIENT_COUNT, (
        f"{n_patients} unique patients among ROIs with labeled cells, "
        f"expected {EXPECTED_PATIENT_COUNT}"
    )
    n_tumor_patients = clinical_with_cells.loc[
        clinical_with_cells["is_tumor"] == "yes", "pat_id"
    ].nunique()
    assert n_tumor_patients == EXPECTED_TUMOR_PATIENT_COUNT, (
        f"{n_tumor_patients} unique is_tumor=='yes' patients, "
        f"expected {EXPECTED_TUMOR_PATIENT_COUNT}"
    )

    # ROI counts do NOT fully reconcile even after the restriction above: 534
    # ROIs have labeled cells (paper: 523 "high-quality" ROIs) and 476 of
    # those are is_tumor=='yes' (paper: 459 "tumor-containing" ROIs). Logged,
    # not asserted — see REPRODUCIBILITY.md Known discrepancies.
    n_rois = len(clinical_with_cells)
    n_tumor_rois = (clinical_with_cells["is_tumor"] == "yes").sum()
    logger.info(f"ROIs with labeled cells: {n_rois} (paper: 523); of those is_tumor=='yes': {n_tumor_rois} (paper: 459)")

    # %% normalized intensities (arcsinh + 99.9th pct censor + min-max, same as clustering input)
    intensity_normalized = prepare_data(base_dir=base_dir, mask_version="annotated")
    assert len(intensity_normalized) == len(metadata), (
        f"intensity rows ({len(intensity_normalized)}) != metadata rows ({len(metadata)})"
    )

    # %%
    metadata.to_parquet(export_dir / "metadata.parquet")
    clinical.to_parquet(export_dir / "clinical.parquet")
    intensity_normalized.to_parquet(export_dir / "intensity_normalized.parquet")
    logger.info(f"Exported metadata/clinical/intensity_normalized to {export_dir}")


if __name__ == "__main__":
    CLI(main)
