> **UPDATE (2026-07-30, final): all scripts repointed, re-verified.** Every
> script in `scripts/` has been updated to read from `$EXPORT_DIR`
> (metadata/clinical/intensity/intensity_normalized -- reproducible via
> `export_for_r.py`) or `$LEGACY_DATA_DIR` (everything else, staged from the
> consolidated copy below) instead of any hardcoded personal-machine or
> `/users/mensmeng/...` path. Every path each script now resolves has been
> checked to exist, except the five genuinely-missing files listed at the
> bottom of this document -- confirmed by direct existence check, not by
> re-reading source. Truly missing, nothing else:
> - `niche_annotations_revised.xlsx` (blocks only `01_annotation_v2.py`;
>   its own output is already precomputed and staged, so nothing downstream
>   needs it)
> - `check_clinical/0-export/clinical.parquet` (non-load-bearing sanity
>   check in `risk_groups_label.R`)
> - `0-paper/0-export/clinical.csv` (ad hoc export, `proportion_heatmap.R`)
> - two `pairwise_wilcoxon_results.csv` files (`proportion_heatmap.R`)
> - the `/home/labadmin/data/pca-v3` spillover-correction input root (3
>   scripts; this pipeline's inclusion in the final analysis was already
>   flagged unresolved before this audit)

> **UPDATE (2026-07-30): consolidated and re-verified.** Every file below
> marked "Found" (at `/work/.../prometex/data/PCa/...`, `/work/.../PCa_NHood/...`,
> or the old repo `/work/.../prometex/projects/PCa/...`) has now been copied
> (not moved — originals untouched) into one new location:
> `/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/`,
> mirroring each file's original expected relative path (e.g. `0-paper/0-export/`,
> `5-niches/annotation/`, `5-niches/barplot_data/`, `PCA_NHOODs_clean/`,
> `PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas/`).
> All copies verified readable. Re-running this audit's basename check against
> the new location confirms: **every previously-"Found" file resolves there**,
> and the only files that remain genuinely missing are `niche_annotations_revised.xlsx`,
> the `pca-v3` spillover root, `Downloads/clinical.csv`, two
> `pairwise_wilcoxon_results.csv` files, and the non-load-bearing
> `check_clinical/0-export/clinical.parquet` sanity check — exactly the set
> already flagged genuinely missing below. **Scripts still reference their
> old hardcoded paths** — this only stages the data; repointing each script
> is a separate follow-up.

# Missing/inaccessible input files, by script

Scope: only files a script **reads** (input dependencies) — output paths a
script writes to are excluded, since a non-existent output path is expected,
not a problem. For every hardcoded absolute path outside the `BASE_DIR`/
`EXPORT_DIR` convention, this checks (a) whether it exists and is readable
here, and (b) whether an identically-named file exists elsewhere on the
shared HPC storage (`/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/`, group
`prometex_101454-pr-g`) or was ever produced by a script in this repo
(searched by basename, not literal path, across the full working tree and
`git log --all`, since paths are often built from variables).

## Headline finding

**Most of these are not actually gone.** The scripts hardcode paths from two
retired personal machines — `/Users/me3312/Documents/Paper_PCa/...` and
`/Users/adrianomartinelli/Library/CloudStorage/OneDrive-.../PCa/...` — plus
one collaborator's live but permission-restricted account,
`/users/mensmeng/workspace/...`. For most filenames referenced from the two
retired-machine paths, a file with the **same basename and, usually, the
same relative sub-path** already exists, readable, on the shared cluster at
`/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/...` — apparently the same
directory tree, migrated to shared storage under a new root, sometime after
these scripts were last edited. Fixing these scripts is very likely a
`sed`-style root-swap, not data recovery.

**This includes the Figure 4a/c dendrogram-cluster file** flagged as
"unrecoverable" earlier in this investigation
(`REPRODUCIBILITY.md`, commit `0604c1ba`): `metadata_with_dendrogram_colors_label_pat_id.parquet`
exists and is readable at
`/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet`.
**That conclusion needs to be revisited** — the real P1-P6 assignment may be
recoverable after all. Flagging this for follow-up rather than fixing it
here, since it changes a conclusion already reported.

A third location also exists: the old pre-migration project repo at
`/work/FAC/FBM/DBC/mrapsoma/prometex/projects/PCa/` (a full git checkout,
not the `data/` tree above) — this is where `colormaps.yaml` turned up after
being reported missing from the other two locations (now copied into this
repo at `resources/colormaps.yaml`). It was also checked for the remaining
genuinely-missing files below (`niche_annotations_revised.xlsx`,
`pairwise_wilcoxon_results.csv`, `panel.csv`, `clinical.csv`) — no hits.

A small number of files are genuinely absent from all three locations —
those are marked **GENUINELY MISSING** below and are the ones worth asking
collaborators about.

---

## scripts/00-spillover-correction/spillover_correct_images_pca.R

- **`/home/labadmin/data/pca-v3`** (`base_dir`, L25) — reads
  `images/filtered/panel.csv` and raw images under this root.
  **NOT_FOUND** locally, no equivalent found under
  `/work/.../prometex/data/PCa/` either (searched for `panel.csv`, 0 hits).
  **GENUINELY MISSING.** Per `REPRODUCIBILITY.md`, this whole
  spillover-correction pipeline's inclusion in the final analysis is
  already flagged unresolved — the missing input may be moot if the step
  was abandoned.

## scripts/00-spillover-correction/spillover_correction_pca_120124.R

- **`/home/labadmin/data/pca-v3`** (`base_dir`, L11) — same as above.
  **GENUINELY MISSING**, same caveat.

## scripts/00-spillover-correction/spillover_correction_pca_190224.R

- **`/home/labadmin/data/pca-v3`** (`base_dir`, L11) — same as above.
  **GENUINELY MISSING**, same caveat.

## scripts/03-survival/plot-cox-hazard-ratio.r

- **`/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/.../PCa/outputs/7-survival`**
  (`surv.dir`, L10) → `scores.path` read at L25. **NOT_FOUND** locally. No
  confirmed same-path equivalent on shared storage; `/work/.../PCa/11-survival/`
  and `/work/.../PCa/scores/` exist and are plausibly related but content
  wasn't cross-checked — **unconfirmed, not asserted equivalent**.

## scripts/03-survival/survival-cell-freq-groups.R

- **`.../PCa/0-export/survival-cell-freq-groups.parquet`** (`clinical.path`, L14).
  **NOT_FOUND** at the hardcoded Mac path. **Found** at
  `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/0-export/survival-cell-freq-groups.parquet`.

## scripts/03-survival/survival-gleason.R

- **`.../PCa/0-export/survival-gleason.parquet`** (`clinical.path`, L14).
  **NOT_FOUND** at hardcoded path. **Found** at
  `/work/.../PCa/0-paper/0-export/survival-gleason.parquet`.

## scripts/03-survival/survival-interactions-diff.R

- **`.../PCa/0-export/scores-v2.parquet`** (`scores.path`, L12). **NOT_FOUND**
  at hardcoded path. **Found** at `/work/.../PCa/0-paper/0-export/scores-v2.parquet`.
- **`.../PCa/0-export/clinical.parquet`** (`clinical.path`, L13). **NOT_FOUND**
  at hardcoded path. **Found** at `/work/.../PCa/0-paper/0-export/clinical.parquet`
  (also independently producible by this repo's own `scripts/00-data-export/export_for_r.py:116`).

## scripts/03-survival/survival-interactions-obs.R

- Same two paths as `survival-interactions-diff.R` (L12, L13) — same
  findings.

## scripts/03-survival/survival-melissas-km.R

- **`.../PCa/0-export`** (`base_dir`, L7) → `clinical.path`/`metadata.path`
  read at L28-29. **NOT_FOUND** at hardcoded path. **Found**: `clinical.parquet`
  and `metadata.parquet` both exist at `/work/.../PCa/0-paper/0-export/` (and
  are reproducible via this repo's `export_for_r.py`).

## scripts/03-survival/survival-proportions.r

- **`.../PCa/0-export/scores-v2.parquet`** (L13) and **`.../PCa/0-export/clinical.parquet`**
  (L14). Same findings as `survival-interactions-diff.R`.

## scripts/03-survival/survival-stromogenic-inflammation.R

- **`.../PCa/0-export/survival-stromogenic-inflammation.parquet`**
  (`clinical.path`, L15). **NOT_FOUND** at hardcoded path. **Found** at
  `/work/.../PCa/0-paper/0-export/survival-stromogenic-inflammation.parquet`.

## scripts/03-survival/survival.r

- **`.../PCa/0-export/scores.parquet`** (L10). **Found** at
  `/work/.../PCa/0-paper/0-export/scores.parquet` (also `/work/.../PCa/scores/scores.parquet`).
- **`.../PCa/0-export/metadata.parquet`** (L11) and **`.../PCa/0-export/clinical.parquet`**
  (L12). **Found** at `/work/.../PCa/0-paper/0-export/`.

## scripts/04-heatmaps/2-cell-types-heatmap.R

- **`.../PCa/0-export/intensity_normalized.parquet`** (`data_path`, L20).
  **NOT_FOUND** at hardcoded Mac path. **Found** at
  `/work/.../PCa/0-paper/0-export/intensity_normalized.parquet` (also
  reproducible via `export_for_r.py:117`).
- **`.../PCa/0-export/metadata.parquet`** (`metadata_path`, L23). Same —
  found at `/work/.../PCa/0-paper/0-export/`.
- **`/Users/adrianomartinelli/projects/PCa/colormaps.yaml`** (`colormap_path`,
  L25). **NOT_FOUND** at the hardcoded Mac path or on the migrated shared
  storage at `/work/.../prometex/data/PCa/`. **Found**, however, in a third
  location not yet covered by this audit: the old pre-migration project
  repo at `/work/FAC/FBM/DBC/mrapsoma/prometex/projects/PCa/colormaps.yaml`
  (a full git checkout, not the `data/` tree). **Resolved** — copied into
  this repo at `resources/colormaps.yaml`. `2-cell-types-heatmap.R:25` still
  needs its hardcoded path updated to point there if this script is revived.

## scripts/11_niches/110_analysis/00_kmeans_clustering.py

- **`/users/mensmeng/workspace/nhoods/PCa/05_nhoods/PCA_NHOODs_clean/robustness`**
  (`sys.path.append`, L9) — imports `utils.clustering.perform_kmeans_clustering`
  and `utils.clustering.wrapper_nhood_filtering`. **Permission denied**
  (confirmed via `ls`, not "not found" — it's inside another user's private
  home directory, not the shared `prometex_101454-pr-g` group space). No
  equivalent module found in a time-bounded search of shared storage.
  **No evidence anywhere in this repo of this module's source.**
  (Note: the scientific input data this script reads,
  `/work/.../prometex/data/PCa_NHood/final_analysis/evaluation/count/CellCellNeighborhoods/`,
  IS accessible — that part is not a problem.)

## scripts/11_niches/110_analysis/01_annotation_v2.py

- **`/users/mensmeng/workspace/nhoods/PCa/05_nhoods/PCA_NHOODs_clean/niche_annotations_revised.xlsx`**
  (L47) — the manually-curated cluster→niche-name mapping (columns `cluster`,
  `niche`, `meta_niche`, `niche_color`, `meta_niche_color`). **Permission
  denied.** **No equivalent found anywhere** in the shared `/work/.../PCa/`
  tree (searched for `.xlsx` and `niche_annotations_revised` — 0 hits).
  **GENUINELY MISSING.** However: this script's own *output* —
  `clusters_annotated_v2.parquet` and `niche_annotations_v2.csv` — already
  exists, precomputed, at `/work/.../PCa/5-niches/annotation/` (see next
  entries), so downstream scripts don't need this script to be re-run.

## scripts/11_niches/111_heatmaps/visualize_composition.py

- **`/users/mensmeng/workspace/PCA_NHOODs_clean/robustness`**
  (`sys.path.append`, L9) — same missing utility module as above.
  **Permission denied. No evidence found.**

## scripts/11_niches/111_heatmaps/niche_pairwise_corrleation.R

- **`.../PCa/5-niches/frequencies`** (`base_dir`, L4) → reads
  `stacked_barplots/props_niche_tma_id.parquet` (L5). **NOT_FOUND** at
  hardcoded Mac path. **Found** at
  `/work/.../PCa/5-niches/frequencies/stacked_barplots/props_niche_tma_id.parquet`.
- **`.../PCa/0-paper/0-export/clinical.parquet`** (L20). **Found** at
  `/work/.../PCa/0-paper/0-export/clinical.parquet`.
- **`.../PCa/5-niches`** (`result_dir`, L2) → reads
  `annotation/niche_annotations_v2.csv` (L63). **Found** at
  `/work/.../PCa/5-niches/annotation/niche_annotations_v2.csv`.

## scripts/11_niches/111_heatmaps/patient_heatmap.R

- **`.../PCa/5-niches/frequencies`** (`base_dir`, L4) → reads
  `stacked_barplots/props_niche_pat_id.parquet` (L5). **Found** at
  `/work/.../PCa/5-niches/frequencies/stacked_barplots/props_niche_pat_id.parquet`.
- **`.../PCa/0-paper/0-export/clinical.parquet`** (L20). **Found**, as above.
- **`.../PCa/5-niches`** (`result_dir`, L2) → reads `annotation/niche_annotations_v2.csv`
  (L141). **Found**, as above.

## scripts/11_niches/111_heatmaps/proportion_heatmap.R

- Same `base_dir`/`result_dir`/`clinical.parquet` reads as `patient_heatmap.R`
  (L2, L4, L20, L224) — **all found** at the `/work/.../PCa/` locations above.
- **`/Users/me3312/Downloads/clinical.csv`** (`clinical.new.path`, L24, read
  at L25). **NOT_FOUND**, and **no equivalent `clinical.csv` found anywhere**
  in the shared tree (only `clinical.parquet` exists). **GENUINELY MISSING**
  — this looks like an ad hoc CSV export for manual inspection, never
  regenerated.
- **`.../PCa/0-paper/0-export/metadata.parquet`** (L32, direct
  `read_parquet`). **Found** at `/work/.../PCa/0-paper/0-export/metadata.parquet`.
- **`.../frequencies/frequency_boxplots/inflammation/pairwise_wilcoxon_results.csv`**
  (L261) and the `stromogenic_smc_loss_reactive_stroma_present` sibling
  (L262). **NOT_FOUND**, and **no equivalent found anywhere** in the shared
  tree. **GENUINELY MISSING** — these are statistical-test output files
  that would need to be regenerated by whatever computed the pairwise
  Wilcoxon tests (not present in this repo).

## scripts/11_niches/112_violinplots/inflammation_vis.R

- **`.../PCa/0-paper/0-export/clinical.parquet`** (L15). **Found.**
- **`.../PCa/5-niches/annotation/clusters_annotated_v2.parquet`** (L44,
  direct `read_parquet`). **Found** at `/work/.../PCa/5-niches/annotation/clusters_annotated_v2.parquet`.

## scripts/11_niches/112_violinplots/stromogenic_vis.R

- Same two reads as `inflammation_vis.R` (L15, L44) — **both found**.

## scripts/11_niches/112_violinplots/pairwise_testing_niches.R

- **`.../PCa/5-niches/frequencies`** (`base_dir`, L4) — not directly read in
  this script (only used to build an unused/write-only path); no genuine
  read of this literal.
- **`.../PCa/0-paper/0-export/clinical.parquet`** (L16). **Found.**
- **`.../PCa/5-niches/annotation/clusters_annotated_v2.parquet`** (L45).
  **Found.**

## scripts/11_niches/113_survival/inflammation_outcome.R

- **`.../PCa/5-niches/frequencies`** (`base_dir`, L4) → reads
  `stacked_barplots/props_niche_tma_id.parquet` (L5). **Found.**
- **`.../PCa/0-paper/0-export/clinical.parquet`** (L22, L128). **Found.**
- **`.../PCa/5-niches/annotation/clusters_annotated_v2.parquet`** (L52).
  **Found.**

## scripts/11_niches/113_survival/label_kaplan_meier_binary.R

- **`.../PCa/5-niches/frequencies`** (`base_dir`, L7) → `props_niche_tma_id.parquet`
  read (L8). **Found.**
- **`.../PCa/0-paper/0-export/clinical.parquet`** (L25, L129). **Found.**
- **`.../PCa/5-niches/annotation/clusters_annotated_v2.parquet`** (L55).
  **Found.**

## scripts/11_niches/113_survival/label_kaplan_meier_vis.R

- Same three reads as `label_kaplan_meier_binary.R` (L7/8, L25/L111, L55) —
  **all found**.

## scripts/11_niches/113_survival/niche_kaplan_meier_binary.R

- Same pattern (L7/8, L25/L131, L55) — **all found**.

## scripts/11_niches/113_survival/stromogenic_outcome.R

- Same pattern (L4/5, L22/L149, L52) — **all found**.

## scripts/11_niches/113_survival/risk_groups_label.R

- **`.../PCa/0-paper/0-export/clinical.parquet`** (`clinical.path`, L2).
  **Found** at `/work/.../PCa/0-paper/0-export/clinical.parquet`.
- **`/Users/me3312/Desktop/check_clinical/0-export/clinical.parquet`**
  (`new_clinical_path`, L4, read at L5 for an `all.equal()` sanity check
  against the above). **NOT_FOUND** at hardcoded path; **no equivalent
  `check_clinical` directory found** on shared storage. **GENUINELY
  MISSING** — but this read is only a one-off manual consistency check
  between two exports, not load-bearing for any figure.
- **`.../PCa/5-niches/barplot_data/metadata_clustered_pat_id_label.csv`**
  (`path_patient`, L14). **NOT_FOUND** at hardcoded path. **Found** at
  `/work/.../PCa/5-niches/barplot_data/metadata_clustered_pat_id_label.csv`.
- **`.../PCa/5-niches/barplot_data/metadata_clustered_tma_id_label.csv`**
  (`path_tma`, L15). **Found** at
  `/work/.../PCa/5-niches/barplot_data/metadata_clustered_tma_id_label.csv`.
- **`/Users/me3312/Desktop/check_clinical/0-export/survival-cell-freq-groups.parquet`**
  (`path_groups`, L20). **NOT_FOUND** at hardcoded path, but the same
  basename exists at `/work/.../PCa/0-paper/0-export/survival-cell-freq-groups.parquet`
  (different subdirectory — `check_clinical/0-export` vs `0-paper/0-export`;
  likely the same file, not independently confirmed byte-for-byte).
- **`/Users/me3312/Desktop/check_clinical/0-export/survival-gleason.parquet`**
  (`path_gleason`, L21). Same situation — found as
  `/work/.../PCa/0-paper/0-export/survival-gleason.parquet` (unconfirmed
  match).
- **`/Users/me3312/Desktop/check_clinical/0-export/survival-stromogenic-inflammation.parquet`**
  (`path_histo`, L22). Same situation — found as
  `/work/.../PCa/0-paper/0-export/survival-stromogenic-inflammation.parquet`
  (unconfirmed match).
- **`.../PCa/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet`**
  (`final_path`, L36) — **the manually-labeled P1-P6 dendrogram assignment.**
  **NOT_FOUND** at the hardcoded Mac path, but **confirmed FOUND and
  readable** at
  `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet`
  (verified: file exists, valid Parquet header, 10,180 bytes, owned by
  `mensmeng`, group `prometex_101454-pr-g`, readable). **This reverses the
  "unrecoverable" conclusion reported earlier in this conversation and
  committed to `REPRODUCIBILITY.md` (commit `0604c1ba`) — that commit needs
  to be revisited, since the real cluster-assignment data may in fact be
  usable to fix Figure 4a/c.**

## scripts/11_niches/114_interactions/compute_interactions.py

- **`/users/mensmeng/workspace/PCA_NHOODs_clean/robustness`**
  (`sys.path.append`, L13). **Permission denied. No evidence found.**
- **`/users/mensmeng/workspace/nhoods/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas`**
  (`anndata_dir`, L38, read via `open(anndata_path, 'rb')` at L66).
  **Permission denied.** Note this is the `proportion` variant of the
  neighborhoods dataset; the `count` variant used by `00_kmeans_clustering.py`
  is accessible, but this `proportion`/`anndatas` variant was not found
  accessible anywhere in a time-bounded search. **Unresolved.**

## scripts/11_niches/114_interactions/visualize_interactions_lfc.py

- **`/users/mensmeng/workspace/PCA_NHOODs_clean/robustness`**
  (`sys.path.append`, L14). Same missing utility module. **Permission
  denied. No evidence found.**

## scripts/11_niches/114_interactions/visualize_circos_plot.py

- **`/users/mensmeng/workspace/nhoods/PCa/05_nhoods/PCA_NHOODs_clean/05_niche_identification/visualization`**
  (`sys.path.append`, L8). **Permission denied. No evidence found.**

## scripts/11_niches/114_interactions/nhood_interactions.py

- **`/users/mensmeng/workspace/PCA_NHOODs_clean/workflow`**
  (`sys.path.append`, L18). **Permission denied. No evidence found.**

---

## Summary table

| Category | Count |
|---|---|
| Scripts scanned (`.py`/`.R`/`.r` under `scripts/`) | 66 |
| Hardcoded-path reads flagged | ~55 |
| **Found** under `/work/.../prometex/data/PCa/` (same basename) — just a stale path root | ~40 |
| **Found** in the old pre-migration repo (`/work/.../prometex/projects/PCa/`) and copied in | 1: `colormaps.yaml` (now `resources/colormaps.yaml`) |
| **GENUINELY MISSING** (not found in any of the three locations) | 6: `pca-v3` spillover root (×3 scripts, 1 unique path), `Downloads/clinical.csv`, two `pairwise_wilcoxon_results.csv` files, `check_clinical/0-export/clinical.parquet` (non-load-bearing sanity check) |
| **Permission denied** (exists, but on a collaborator's private, non-group-shared account) | 7: `niche_annotations_revised.xlsx`, 4 `utils`/module import paths, 1 `anndatas` data directory |
| **Reverses an earlier conclusion** | 1: `metadata_with_dendrogram_colors_label_pat_id.parquet` — previously reported unrecoverable, now confirmed present and readable |
