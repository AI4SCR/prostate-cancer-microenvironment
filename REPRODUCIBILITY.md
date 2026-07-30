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
cp .env.example .env   # fill in BASE_DIR and EXPORT_DIR
uv sync
```

`BASE_DIR` (the `ai4bmr_datasets.PCa` staging root) is **read-only** — no
script in this repo may ever write there, and `assert_outside_base_dir()` /
`resolve_export_dir()` in `src/prostate_cancer/utils.py` enforce that at
runtime (checked via `os.path.realpath`, so a symlinked alias into the same
tree is caught too — not just the literal `BASE_DIR` path). `EXPORT_DIR` is
where every generated table and figure goes instead; point it somewhere
inside this repo (e.g. `data/0-export`, the `.env.example` default) so
everything produced by these scripts lives together and is easy to inspect
or delete.

R scripts load the same `.env` via `dotenv::load_dot_env()` and
`Sys.getenv("BASE_DIR")` / `Sys.getenv("EXPORT_DIR")`. Required R packages:
`dotenv`, `arrow`, `tidyverse`, `survival`, `ggsurvfit`, `gtsummary`,
`compositions`, `coxme`, `ComplexHeatmap`. No `renv.lock` is maintained —
package versions aren't pinned on the R side.

`ai4bmr-datasets` is pinned in `pyproject.toml` to its `pca` git branch: the
`PCa` dataset class this repo depends on currently only exists there, not on
`ai4bmr-datasets`'s `main` branch (it was removed from `main` in a later
commit). Revisit the pin once/if `pca` merges upstream.

## Data layout

`ai4bmr_datasets.PCa` stages the read-only dataset under `$BASE_DIR`:

```
$BASE_DIR/
├── 01_raw/
│   ├── raw/                      # input .mcd files (external)
│   ├── acquisitions/             # process_acquisitions() output: per-acquisition TIFF+JSON
│   ├── masks/deepcell/           # EXTERNAL: segmentation masks must already be staged here
│   ├── clustering/annotations.parquet   # produced by scripts/01-clustering/09_main-annotate.py
│   ├── reclustering/, reclustering-v2/  # manual re-clustering memberships (external/manual)
│   └── metadata/                 # label-names.xlsx, ROI_matching_blockID.xlsx, tma-annotations-v3.xlsx
└── 02_processed/
    ├── images/{raw,filtered}/
    ├── masks/{deepcell,filtered,annotated}/
    ├── metadata/{clinical.parquet, filtered-annotated/*.parquet}
    └── features/{intensity,spatial}/{image_version}-{mask_version}/*.parquet
```

Everything this repo generates goes to `$EXPORT_DIR` instead (outside
`$BASE_DIR`):

```
$EXPORT_DIR/
├── metadata.parquet, clinical.parquet, intensity_normalized.parquet   # export_for_r.py
└── figures/figureN/*.png                                             # scripts/figures/figureN_*.py
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
   `metadata.parquet`, `clinical.parquet`, `intensity.parquet`,
   `intensity_normalized.parquet` to `$EXPORT_DIR` for the R scripts (see
   below). This is a 1:1 port of the original publication's own export
   script (see "Normalization: `prepare_data()` vs `normalize()`" below) —
   verified byte-identical against the legacy exports.
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
| Fig 4 | (a) patient composition clusters P1–P6 (JSD hierarchical clustering, max-pooled proportions), (b) metagroup distribution, (c) KM by patient group, (d-e) Cox PH (overall survival / progression) | Final annotated cell table | `scripts/figures/figure4_patient_clustering.py` (a, reconstructed — see Known discrepancies for a confirmed mismatch on the resulting clusters), `figure4_survival.R` (c-e — d/e match the paper exactly, c does not) |
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
- **ROI / patient count — RESOLVED.** The paper reports 523 "high-quality"
  ROIs (Results) and, after further removing ROIs lacking clinical
  annotation or tumor content, 459 "tumor-containing" ROIs (Methods, "IMC
  Data Acquisition") — two sequential QC stages, not a contradiction.

  **Investigated 2026-07-29.** `sample_id` (the ai4bmr-datasets ROI key) is
  **per acquisition**, not per physical core: some cores had an interrupted
  IMC scan that was re-acquired, producing two `sample_id` rows for the same
  physical ROI (e.g. `231209_ibl_x3y19_198_1_13` and
  `231209_ibl_x3y19_198_12`, both `pat_id` 96.3799, same `tma_id`
  `IBL_X3Y19` — one annotated `notes = "lot missing"`). `tma_id` is the
  physical-core identifier the paper actually counts. Deduplicating by
  `tma_id` reconciles every number **exactly**:
  - All 542 clinical rows → **523** unique `tma_id` = paper's "523
    high-quality ROIs", exactly.
  - The 476 rows with labeled cells AND `is_tumor == "yes"` → **459** unique
    `tma_id` = paper's "459 tumor-containing ROIs", exactly.
  - (The 534 rows with labeled cells, pre-tumor-filter → 515 unique
    `tma_id` — this number was already known but undocumented: it's a stray
    comment, `# 515 TMA_IDs before filtering for tumor cores only`, in
    `scripts/03-survival/survival-proportions.r`.)

  Patient counts (previous paragraph) are also explained by the same
  `sample_id`-vs-`tma_id` distinction plus one ROI that was never processed
  into cells.

  **Implication for every figure branch**: any ROI/core-level aggregation
  (cell counts per core, niche abundance per core, KM stratification by
  core) must group by `tma_id`, not `sample_id`, or interrupted-acquisition
  duplicates will be double-counted. `export_for_r.py`'s exported
  `clinical.parquet` is restricted to ROIs present in both `metadata` and
  `clinical` (534 rows, matching the original publication export script) —
  so it asserts **515** unique `tma_id` (the "534 rows, pre-tumor-filter"
  case above), not the paper's headline 523, which describes the full
  unrestricted 542-row clinical table this repo never exports on its own.
  459 (tumor-containing ROIs among labeled cells) is still asserted as before.

  **Authoritative reference**, confirming 515/523/459 are three different
  cohort splits, not competing values for the same thing — from the old
  repo's own `000_paper/0-export/readme.md`:

  | Split | #patients | #acquisitions | #ROIs | #tumor | #no_tumor | #N/A |
  |---|:-:|:-:|:-:|:-:|:-:|:-:|
  | Clinical (full, unrestricted) | 196 | 542 | **523** | 481 | 40 | 21 |
  | Clinical annotated (restricted to sample_ids with labeled cells) | 195 | 534 | **515** | 476 | 39 | 19 |
  | Clinical tumor-only | 190 | 476 | **459** | 476 | 0 | 0 |

  `export_for_r.py`'s exported `clinical.parquet` is exactly the "Clinical
  annotated" row (534 rows / 515 ROIs / 195 patients) — hence
  `EXPECTED_ROI_COUNT = 515`, `EXPECTED_PATIENT_COUNT = 195`. The paper's
  headline 523/196 describes the row above it (unrestricted); 459/190 the
  row below (tumor-only, computed separately from `clinical_with_cells` in
  this script). All three are internally consistent — filtering down through
  the same three stages the paper itself describes (patient-level: 196 → 195
  → 190; acquisition/ROI-level: 542→523 → 534→515 → 476→459, each split
  losing 19 ROIs to interrupted-and-restarted acquisitions).
- **Normalization: `prepare_data()` vs `normalize()` — RESOLVED, was a
  wrong-function port, not a data or formula discrepancy.**
  `export_for_r.py` used to build `intensity_normalized.parquet` by calling
  `prepare_data()`. That produced values differing from the legacy
  publication export on ~79% of cells (max abs diff 0.063 on a min-max-scaled
  [0,1] range) -- small, but real and systematic, not floating-point noise.

  **Investigated 2026-07-30.** Traced to the pre-migration repo's own
  `000_paper/0-export/data.py` (the actual script that produced the
  published `0-paper/0-export/` tables): it loads `PCa(mask_version=
  'annotated', ...)` -- same as this repo -- then calls the *separate*
  `utils.normalize(intensity, exclude_zeros=True)`, never `prepare_data()`.
  The old repo had three near-duplicate normalization implementations
  (`normalize_data()`, `normalize()`, and `prepare_data()`'s own inline
  copy); `prepare_data()`'s inline copy was itself dead code in that repo,
  never used to produce any published output, and never excludes zeros when
  computing the 99.9th-percentile censoring threshold. This repo's
  `prepare_data()` was a faithful, unmodified copy of that dead code,
  carried over since this repo's very first commit -- a name-collision bug,
  not a strategy that changed over time (confirmed: `git log -p` on the old
  repo's `utils.py` shows `prepare_data()`'s normalization block and the
  `mask_version="cleaned"`→`"annotated"` question were red herrings; neither
  the formula nor these dataset parameters ever changed there).

  **Fix**: `export_for_r.py` now ports `000_paper/0-export/data.py` itself
  (paths only changed) and calls `utils.normalize(intensity,
  exclude_zeros=True)` directly, matching `mask_version='annotated'`,
  `load_intensity=True`, `align=False`, and the `sample_ids = set(metadata...)
  & set(clinical...)` restriction exactly. **Verified**: regenerated
  `clinical.parquet` (534 rows), `metadata.parquet` (2,191,967 rows),
  `intensity.parquet`, and `intensity_normalized.parquet` (both
  2,191,967×40) are now byte-identical to the legacy exports -- 0 cells
  differing out of 87,678,680, at `atol=1e-9`. `prepare_data()` remains in
  `src/prostate_cancer/utils.py` only for the 01-clustering bootstrap
  scripts (`mask_version="filtered"`), which is a different, still-valid use
  case -- it is no longer called anywhere in the export path.
- **`scripts/01-clustering/03_epithelial-non-epithelial-annotate.py`** is an
  empty file (0 bytes) in the current repo. Not reconstructed here — flagged
  for the `figure-2-cell-phenotyping` branch to investigate.
- **Figure 4a/c patient-composition clustering — UNRESOLVED, confirmed
  mismatch.** `figure4_patient_clustering.py` reconstructs the paper's
  described algorithm exactly as written (per-`tma_id` proportions with a
  pseudocount of 1, max-pooled to one composition vector per patient,
  average-linkage hierarchical clustering on Jensen-Shannon divergence, cut
  to 6 clusters). Run against the real dataset it produces 195 patients
  across P1-P6 with sizes `3 / 101 / 23 / 2 / 65 / 1` — three of the six
  clusters have 1-3 patients. The paper's Fig 4c states Kaplan-Meier found
  "no statistically significant association" between these groups and
  survival; the reconstructed clusters give log-rank p=6.65e-09, driven
  almost entirely by the three tiny clusters' wild early-event swings (a
  3-patient and a 1-patient cluster each show a near-vertical survival drop).
  This is very unlikely to be the same partition the paper reports.
  **Not resolved, and now confirmed unrecoverable from code.** The Methods
  text doesn't specify how the dendrogram was cut into exactly 6 groups (a
  fixed height threshold would very plausibly give a more balanced partition
  than forcing `k=6` via `scipy`'s `maxclust`, which was the choice made
  here) — but this is not just an underspecified Methods section. The only
  other script in this repo that touches a patient-level cluster/group
  assignment, `scripts/11_niches/113_survival/risk_groups_label.R`, reads a
  manually-produced `leaf_color_group` column (`risk_groups_label.R:66`)
  from `metadata_with_dendrogram_colors_label_pat_id.parquet`
  (`risk_groups_label.R:36`) and a `cluster_group` column from
  `metadata_clustered_pat_id_label.csv` (`risk_groups_label.R:14,17,29`) —
  both under a hardcoded local path,
  `/Users/me3312/Documents/Paper_PCa/5-niches/barplot_data/`, that exists on
  neither this repo, its full `git log --all` history (the script entered in
  the initial `init` commit already referencing these paths), nor the
  published Zenodo archive. No script anywhere in this repo's history writes
  either file. The "dendrogram colors label" naming strongly suggests the
  real cluster assignment was produced by manually coloring/labeling
  dendrogram leaves (by eye, likely in a GUI or notebook), not by any single
  reproducible height/k threshold — so no choice of `scipy` cut parameter in
  `figure4_patient_clustering.py` is expected to recover the paper's exact
  P1-P6 partition. This is a genuine gap in what the paper's code release
  preserved, not a bug in this reconstruction.
  **However**: Figure 4d-e (Cox PH per cell type on the same underlying
  CLR-transformed composition vectors, independent of the P1-P6 clustering)
  reproduce the paper's headline finding exactly — `epithelial-luminal(ERG+p53+)`
  is the sole FDR-significant hit for both overall survival and disease
  progression, matching "the luminal ERG+p53+ epithelial population emerged
  as the only cell type significantly associated with both overall and
  progression-free survival." The composition/CLR pipeline itself is
  therefore validated; only the specific dendrogram cut for Fig 4a/c remains
  open.
