# Reproducibility guide

This repository is the code companion to *"Spatial single-cell proteomics defines
multicellular niches in the primary prostate cancer microenvironment"*
(bioRxiv 10.64898/2026.04.30.721907v1). The paper's Code Availability statement
promises that all scripts needed to reproduce the main figures live here. This
document is the map: the full pipeline from raw acquisitions to each figure, in
run order, plus every known gap and discrepancy — documented rather than
silently patched, per this project's conventions (see `CLAUDE.md`).

## Setup

```bash
cp .env.example .env   # fill in BASE_DIR (and EXPORT_DIR if not colocated with BASE_DIR)
uv sync
```

R scripts load the same `.env` via `dotenv::load_dot_env()` and
`Sys.getenv("BASE_DIR")` / `Sys.getenv("EXPORT_DIR")`. Required R packages:
`dotenv`, `arrow`, `tidyverse`, `survival`, `ggsurvfit`, `gtsummary`,
`compositions`, `coxme`, `ComplexHeatmap`. No `renv.lock` is maintained —
package versions aren't pinned on the R side.

`ai4bmr-datasets` is pinned in `pyproject.toml` to its `pca` git branch: the
`PCa` dataset class this repo depends on currently only exists there, not on
`ai4bmr-datasets`'s `main` branch (it was removed from `main` in a later
commit). Revisit the pin once/if `pca` merges upstream.

## Data layout (`ai4bmr_datasets.PCa`, staged under `$BASE_DIR`)

```
$BASE_DIR/
├── 01_raw/
│   ├── raw/                      # input .mcd files (external)
│   ├── acquisitions/             # process_acquisitions() output: per-acquisition TIFF+JSON
│   ├── masks/deepcell/           # EXTERNAL: segmentation masks must already be staged here
│   ├── clustering/annotations.parquet   # produced by scripts/01-clustering/09_main-annotate.py
│   ├── reclustering/, reclustering-v2/  # manual re-clustering memberships (external/manual)
│   └── metadata/                 # label-names.xlsx, ROI_matching_blockID.xlsx, tma-annotations-v3.xlsx
├── 02_processed/
│   ├── images/{raw,filtered}/
│   ├── masks/{deepcell,filtered,annotated}/
│   ├── metadata/{clinical.parquet, filtered-annotated/*.parquet}
│   └── features/{intensity,spatial}/{image_version}-{mask_version}/*.parquet
└── 0-export/                     # EXPORT_DIR: R-facing parquet exports (see below)
```

## Pipeline: raw acquisitions → labeled cells (run once, in order)

Segmentation is **not** performed by `ai4bmr-datasets`. Steinbock's hot-pixel
removal and DeepCell nuclear segmentation (`steinbock preprocess imc images
--hpf 50`, `steinbock segment deepcell --minmax --type nuclear`, matching the
Methods section and the commands already documented in `README.md`) must be
run externally, with the resulting masks staged at `01_raw/masks/deepcell/`
before step 2 below.

1. **External**: raw `.mcd` files under `01_raw/raw/`; deepcell masks under
   `01_raw/masks/deepcell/`.
2. `ds.process_acquisitions()` → `create_panel()` → `create_images()` →
   `create_filtered_images()` → `create_masks()` → `create_filtered_masks()` →
   `create_clinical_metadata()` → `create_tma_annotations()`.
3. `ds.compute_features(image_version="filtered", mask_version="filtered")` —
   features for **all** segmented cells, no cell-type labels yet. This is the
   correct input for clustering on a from-scratch `BASE_DIR` (as opposed to
   `mask_version="annotated"`, which requires labels that don't exist yet —
   see "the bootstrap loop" below).
4. `scripts/01-clustering/01_immune-non-immune.py` through
   `09_main-annotate.py`, in numeric order. Each stage clusters one cellular
   compartment (immune / epithelial / stromal / endothelial / undefined /
   basal) using `src/prostate_cancer/cluster.py:cluster()` and the matching
   `*-annotate.py` script's marker-based manual labeling. All of these now
   call `prepare_data(..., mask_version="filtered")` — see below for why.
5. `09_main-annotate.py` writes the combined manual annotations to
   `01_raw/clustering/annotations.parquet` — the exact path
   `PCa.label_transfer()` reads (this was previously a path mismatch; fixed).
6. `ds.label_transfer()` → `ds.create_annotated()` →
   `ds.compute_features(image_version="filtered", mask_version="annotated")`
   → `ds.create_meta_labels()`. This produces the final per-cell labeled
   feature table (`02_processed/features/{intensity,spatial}/filtered-annotated/`
   and `02_processed/metadata/filtered-annotated/`) that every figure script
   consumes.
7. `python scripts/00-data-export/export_for_r.py` — writes
   `metadata.parquet`, `clinical.parquet`, `intensity_normalized.parquet` to
   `$EXPORT_DIR` for the R scripts (see below).
8. Figure branches (Fig 2–7) consume the outputs of steps 6–7.

### The bootstrap loop, and why `mask_version` matters

`PCa.label_transfer()` reads `01_raw/clustering/annotations.parquet` — a file
this repo's own clustering scripts produce, not `ai4bmr-datasets`. And the
clustering scripts get their input intensities from
`src/prostate_cancer/utils.py:prepare_data()`, which used to hardcode
`mask_version="annotated"`. That only worked because the already-published
dataset already has `labels.parquet` materialized — on a genuinely fresh
`BASE_DIR` it's circular (`"annotated"` doesn't exist until label transfer
runs, and label transfer needs the clustering scripts' output first).

`prepare_data()` now takes `mask_version` as a parameter. Step 4 above
(clustering, bootstrap) uses `mask_version="filtered"`; step 6 onward (label
transfer and everything downstream) uses `mask_version="annotated"` (the
default). If you're re-running clustering against an already-labeled
`BASE_DIR` (e.g. the published Zenodo archive), `"filtered"` still works —
it's just a superset of cells that includes ones later excluded during label
transfer.

`scripts/01-clustering/10_transfer_labels.py` was a standalone script that
duplicated what `PCa.label_transfer()` now does inside `ai4bmr-datasets`; it
has been removed as a second, divergent source of truth.

## Figure → script mapping

| Figure | Panels | Upstream (must run first) | Scripts |
|---|---|---|---|
| Fig 1 | workflow schematic, panel table, representative images | — | Not code-derived (BioRender schematic + raw image crops); no script found for the representative-image panels — **gap**, not reconstructed. |
| Fig 2 | (a) 34-cell-type heatmap, (b) UMAP of all cells | `01-clustering/01_*.py` … `09_main-annotate.py` | `04-heatmaps/2-cell-types-heatmap.R`; `02-umaps/0-umaps.py`, `0-umaps-main-types.py` |
| Fig 3 | (a) CAF UMAP, (b) CAF subcluster heatmap, (c-e) representative ROIs | `01-clustering/04_stroma*.py` | `02-umaps/0-umaps-cafs.py`; panel (b) may reuse `2-cell-types-heatmap.R` filtered to CAF labels — **verify while implementing**; (c-e) representative-ROI crops — **gap** |
| Fig 4 | (a) patient composition clusters P1–P6 (JSD hierarchical clustering, max-pooled proportions), (b) metagroup distribution, (c) KM by patient group, (d-e) Cox PH (overall survival / progression) | Final annotated cell table | (a) **gap**: no script found producing the patient-level clustering — reconstruct from Methods ("Cell Type Proportion Quantification"); `03-survival/survival.r`, `survival-proportions.r`, `plot-cox-hazard-ratio.r`, `survival-cell-freq-groups.R` |
| Fig 5 | (a) 18-niche z-scored heatmap, (b) niche correlation matrix, (c) representative cores | Final annotated cell table | `11_niches/110_analysis/00_kmeans_clustering.py` (k=24, `k-means++`, 50-seed ARI robustness sweep — Suppl. Fig 7 has **no discovered script**, may need writing), `01_annotation_v2.py`; `111_heatmaps/z_score_heatmap.R`, `niche_pairwise_corrleation.R`; (c) representative cores — **gap** |
| Fig 6 | (a) niche abundance heatmap, (b) KM niche 6, (c) composition barplot, (d) niche vs inflammation, (e) immune-niche risk score KM | Fig 5 niche labels | `111_heatmaps/patient_heatmap.R`/`proportion_heatmap.R`, `112_violinplots/inflammation_vis.R`, `113_survival/niche_kaplan_meier_binary.R`, `label_kaplan_meier_binary.R`, `risk_groups_label.R`, `inflammation_outcome.R` |
| Fig 7 | (a) stromal niche vs stromogenic status, (b-c) KM niche 9 vs myCAF, (d-f) circos interaction plots, (g-h) representative cores | Fig 5 niche labels | `112_violinplots/stromogenic_vis.R`, `pairwise_testing_niches.R`, `113_survival/stromogenic_outcome.R`, `114_interactions/compute_interactions.py`, `circos_plots.py`, `visualize_circos_plot.py`, `nhood_interactions.py`, `visualize_interactions_lfc.py`; (g-h) representative cores — **gap** |

Paper-reported baseline numbers to sanity-check against: 2,191,967 cells / 195
patients (190 in the final analytical cohort) / 34 cell types / 18 niches /
patient clusters P1–P6 / niches 1–18.

## Known discrepancies (documented, not silently resolved)

- **Spillover correction**: the paper's Methods state channel spillover was
  negligible and required "no further corrections." `scripts/00-spillover-correction/`
  nonetheless contains active R scripts (`spillover_correct_images_pca.R` and
  two dated variants) that compute and apply spillover compensation via
  CATALYST, plus a Python step that compresses the "compensated" images.
  **Unresolved**: confirm with the authors whether this step was part of the
  pipeline that produced the published results, or an earlier abandoned
  attempt superseded by the "negligible crosstalk" finding.
- **Cell count**: the paper reports 2,191,967 cells consistently (Abstract,
  Results, Methods). `ai4bmr_datasets.PCa.label_transfer()` asserts
  `len(annotations) == 2214046` after merging in the reclustering-v2
  memberships — a different, larger number. The final "annotated" table
  further drops ~3% of cells reported as "unclassified" in Results, which may
  reconcile the gap. **Verified 2026-07-29** against the materialized dataset
  at `$BASE_DIR`: `02_processed/metadata/filtered-annotated/*.parquet` sums
  to exactly **2,191,967** cells across **534** sample files — matches the
  paper exactly. `scripts/00-data-export/export_for_r.py` hard-asserts this
  number against the final exported table.
- **ROI / patient count — patient counts fully reconcile, ROI counts don't.**
  The paper reports 523 "high-quality" ROIs (Results) and, after further
  removing ROIs lacking clinical annotation or tumor content, 459
  "tumor-containing" ROIs (Methods, "IMC Data Acquisition") — two sequential
  QC stages, not a contradiction (corrected from an earlier note here that
  read them as inconsistent).

  **Investigated 2026-07-29.** `02_processed/metadata/clinical.parquet` has
  542 ROI-level rows / 196 unique `pat_id` — more than the paper's numbers.
  Restricting to the 534 ROIs that actually have processed, labeled cells
  (i.e. present in `metadata`) resolves the patient counts **exactly**:
  195 unique patients (matches the paper's initial cohort), and 190 of those
  with `is_tumor == "yes"` (matches the paper's final analytical cohort —
  195 − 190 = 5, matching the Methods statement that "five patients were
  removed due to the absence of tumor content"). The one extra patient
  (`pat_id` `95.8128`) has a single clinical ROI row (`240219_ivb_x5y2_8_4`,
  `is_tumor == "yes"`) with **no** corresponding entry in the labeled-cell
  dataset — i.e. that ROI has clinical annotation but was never
  processed into cells (failed segmentation or dropped before the cell atlas
  was built).

  The **ROI counts** remain unreconciled even after this restriction: 534
  ROIs have labeled cells vs. the paper's 523, and 476 of those are
  `is_tumor == "yes"` vs. the paper's 459 — an 11- and 17-ROI gap
  respectively. The only ROI-level exclusion found anywhere in the current
  scripts is a single hardcoded sample (`"240223_012"`, excluded in both
  `09_main-annotate.py` and `PCa.label_transfer()`), nowhere near enough to
  account for the gap. Likely explanation: the "regions presenting staining
  artifacts or obvious segmentation errors were removed" QC step (Methods,
  "Data Processing") was applied manually/interactively at some point and
  isn't captured by any flag in the current `clinical.parquet` or by any
  script in this repo. Not resolved — `export_for_r.py` logs both ROI counts
  (not asserted) so this stays visible rather than silently passing.
- **`scripts/01-clustering/03_epithelial-non-epithelial-annotate.py`** is an
  empty file (0 bytes) in the current repo. Not reconstructed here — flagged
  for the `figure-2-cell-phenotyping` branch to investigate.
