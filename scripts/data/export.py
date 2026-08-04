# %%
"""Export the common per-cell / per-patient tables that every R script needs.

R has no ai4bmr-datasets binding, so every R script in this repo reads Parquet
files from `EXPORT_DIR` instead. This script produces the tables shared across
figure branches: `metadata.parquet` (per-cell labels), `clinical.parquet`
(per-ROI/patient clinical annotations), `intensity.parquet` (raw marker
intensities), and `intensity_normalized.parquet` (arcsinh + 99.9th-percentile
censor, zeros excluded from the censoring threshold, + min-max). Figure-
specific score tables (e.g. `scores.parquet`, `survival-*.parquet`, niche
abundance tables) are produced by each figure branch's own scripts, not here.

This is a 1:1 port (paths only changed) of the original publication's export
script, `000_paper/0-export/data.py` in the pre-migration repo -- see
REPRODUCIBILITY.md for why: `src/prostate_cancer/utils.py:prepare_data()`
looked like the right helper to call (same name as this repo's own dead-code
duplicate in the old repo) but is NOT what produced the published
`intensity_normalized.parquet` -- that came from `utils.normalize(...,
exclude_zeros=True)`, called directly from the export script, never from
`prepare_data()`. Confirmed byte-identical against the legacy export; see
REPRODUCIBILITY.md Known discrepancies.

Requires the full labeling pipeline to have already run (see
REPRODUCIBILITY.md) so that `01_raw/annotations/labels.parquet` exists.
"""
from pathlib import Path

from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import (
    NON_MARKER_CHANNELS,
    assert_outside_base_dir,
    normalize,
    resolve_base_dir,
    resolve_export_dir,
)

# The counts below are fixed properties of the specific, already-published
# dataset snapshot this pipeline reproduces (Zenodo 10.5281/zenodo.19552665) —
# not a moving target a growing cohort would outrun. Asserted (not logged) so
# a stale/different BASE_DIR fails loudly instead of silently producing
# figures that don't match the paper. Reported 3x consistently in the paper
# (Abstract, Results, Methods). See REPRODUCIBILITY.md Known discrepancies
# for how each was confirmed against the real data, and for a related count
# found in ai4bmr-datasets' own code.
EXPECTED_CELL_COUNT = 2_191_967
EXPECTED_PATIENT_COUNT = 195  # initial TMA cohort
EXPECTED_TUMOR_PATIENT_COUNT = 190  # final analytical cohort (is_tumor == "yes")
EXPECTED_ROI_COUNT = 515  # unique physical cores (tma_id), clinical restricted to
# sample_ids present in both metadata and clinical (matching the original
# 000_paper/0-export/data.py's `sample_ids = ... & ...` restriction) --
# confirmed against the legacy clinical.parquet. The previous 523 counted
# clinical rows before that restriction, a different (and never-produced) universe.
EXPECTED_TUMOR_ROI_COUNT = 459  # unique tma_id, restricted to labeled cells + is_tumor == "yes"


def main(base_dir: Path | None = None, export_dir: Path | None = None):
    from ai4bmr_datasets import PCa

    base_dir = Path(base_dir).expanduser() if base_dir else resolve_base_dir()
    export_dir = assert_outside_base_dir(Path(export_dir).expanduser()) if export_dir else resolve_export_dir()
    export_dir.mkdir(parents=True, exist_ok=True)

    # %% per-cell labels + per-ROI clinical annotations + intensities, exactly
    # as the original 000_paper/0-export/data.py loaded them (paths only changed)
    ds = PCa(
        base_dir=base_dir,
        image_version="filtered",
        mask_version="annotated",
        load_intensity=True,
        load_metadata=True,
        align=False,
    )
    ds.setup()

    clinical = ds.clinical.copy()
    metadata = ds.metadata.copy()
    intensity = ds.intensity.copy()

    # restrict to ROIs present in both tables, exactly as the original script did
    sample_ids = sorted(set(metadata.index.get_level_values("sample_id")) & set(clinical.index))
    clinical = clinical.loc[sample_ids]
    metadata = metadata.loc[sample_ids]
    intensity = intensity.loc[sample_ids]
    assert len(metadata) == len(intensity)
    metadata, intensity = metadata.align(intensity, axis=0, join="inner")

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

    # `sample_id` is per-acquisition, not per physical core: interrupted scans
    # were re-acquired, producing two `sample_id` rows for the same `tma_id`.
    # The paper counts physical cores (`tma_id`), not acquisitions — see
    # REPRODUCIBILITY.md Known discrepancies for how this was confirmed.
    # The 523 count is over ALL clinical rows (acquired cores, pre-QC); the
    # 459 count is restricted to cores with labeled cells AND is_tumor=="yes".
    assert "tma_id" in clinical.columns, "clinical table is missing tma_id"
    n_rois = clinical["tma_id"].nunique()
    assert n_rois == EXPECTED_ROI_COUNT, f"{n_rois} unique tma_id, expected {EXPECTED_ROI_COUNT}"

    n_tumor_rois = clinical_with_cells.loc[
        clinical_with_cells["is_tumor"] == "yes", "tma_id"
    ].nunique()
    assert n_tumor_rois == EXPECTED_TUMOR_ROI_COUNT, (
        f"{n_tumor_rois} unique is_tumor=='yes' tma_id with labeled cells, expected {EXPECTED_TUMOR_ROI_COUNT}"
    )

    # %% normalized intensities: arcsinh + 99.9th-pct censor (zeros excluded
    # from the censoring threshold) + min-max -- exactly what
    # 000_paper/0-export/data.py called, confirmed byte-identical against the
    # legacy export. NOT prepare_data(), which never produced this table.
    intensity_normalized = normalize(intensity, exclude_zeros=True)

    # %%
    metadata.to_parquet(export_dir / "metadata.parquet")
    clinical.to_parquet(export_dir / "clinical.parquet")
    intensity.to_parquet(export_dir / "intensity.parquet")
    intensity_normalized.to_parquet(export_dir / "intensity_normalized.parquet")

    # R has no equivalent of `prostate_cancer.utils.NON_MARKER_CHANNELS` to
    # import, so hand it the same list as a plain-text sidecar (one per line)
    # instead of letting each R figure script hardcode its own copy. This is
    # a fixed constant, not data derived from a dataset, so it's written to
    # resources/ rather than EXPORT_DIR.
    resources_dir = Path(__file__).resolve().parents[2] / "resources"
    (resources_dir / "non_marker_channels.txt").write_text("\n".join(NON_MARKER_CHANNELS) + "\n")

    logger.info(f"Exported metadata/clinical/intensity/intensity_normalized to {export_dir}")


if __name__ == "__main__":
    CLI(main)
