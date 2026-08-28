from pathlib import Path

from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import (
    NON_MARKER_CHANNELS,
    assert_outside_base_dir,
    normalize,
    resolve_base_dir,
    resolve_data_dir,
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
EXPECTED_ROI_COUNT = 515  # unique physical cores (tma_id), clinical restricted to sample_ids present in both metadata and clinical
EXPECTED_TUMOR_ROI_COUNT = 459  # unique tma_id, restricted to labeled cells + is_tumor == "yes"

# ds.clinical carries 57 columns; only these are ever read by a script in
# this repo (grepped across scripts/ + src/, one pass per column name).
# Everything else -- block/slide identifiers (original_block_number,
# donor_block_id, unique_tma_sample_id_1..4, slide_code, tma_sample_id,
# tma_coordinates), napari/annotation bookkeeping (napari_sample_id,
# file_name_napari, annotation_roi_he), free text (description, notes),
# an unused duplicate id (patient_id, vs. the actually-used pat_id), a
# duplicate never referenced (gleason_pattern_tma_core, vs. gleason_grp),
# and two never-referenced fields (cgs_pat_1, cgs_pat_2) -- is dropped so
# clinical.parquet only carries what's actually consumed downstream.
CLINICAL_COLUMNS = [
    "pat_id", "tma_id",
    "age_at_surgery", "psa_at_surgery",
    "last_fu", "cause_of_death", "os_status",
    "psa_progr", "psa_progr_time",
    "clinical_progr", "clinical_progr_time",
    "disease_progr", "disease_progr_time", "recurrence_loc", "recurrence",
    "cgrading_biopsy", "cgs_score", "gs_pat_1", "gs_pat_2",
    "gs_grp", "gleason_score", "gleason_score_sum", "gleason_grp",
    "ct_stage", "pt_stage", "pgs_score", "ln_status",
    "surgical_margin_status", "adj_adt", "adj_radio", "d_amico_risk",
    "stromogenic_smc_loss_reactive_stroma_present", "non_stromogenic_smc_abundant",
    "inflammation", "glandular_atrophy_pin", "cribriform", "is_tumor",
]


def main(base_dir: Path | None = None, data_dir: Path | None = None):
    from ai4bmr_datasets import PCa

    base_dir = Path(base_dir).expanduser() if base_dir else resolve_base_dir()
    data_dir = assert_outside_base_dir(Path(data_dir).expanduser()) if data_dir else resolve_data_dir()
    cells_dir = data_dir / "cells"
    cells_dir.mkdir(parents=True, exist_ok=True)

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
    clinical = clinical[CLINICAL_COLUMNS]
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

    # %% normalized intensities: arcsinh + 99.9th-pct censor (zeros excluded from the censoring threshold) + min-max
    intensity_normalized = normalize(intensity, exclude_zeros=True)

    # %%
    metadata.to_parquet(cells_dir / "metadata.parquet")
    clinical.to_parquet(data_dir / "clinical.parquet")
    intensity.to_parquet(cells_dir / "intensity.parquet")
    intensity_normalized.to_parquet(cells_dir / "intensity_normalized.parquet")

    # R has no equivalent of `prostate_cancer.utils.NON_MARKER_CHANNELS` to
    # import, so hand it the same list as a plain-text sidecar (one per line)
    # instead of letting each R figure script hardcode its own copy. This is
    # a fixed constant, not data derived from a dataset, so it's written to
    # resources/ rather than DATA_DIR.
    resources_dir = Path(__file__).resolve().parents[2] / "resources"
    (resources_dir / "non_marker_channels.txt").write_text("\n".join(NON_MARKER_CHANNELS) + "\n")

    logger.info(f"Exported metadata/intensity/intensity_normalized to {cells_dir}, clinical to {data_dir}")


if __name__ == "__main__":
    CLI(main)
