# Bugs found and fixed while executing the ported figure scripts

Audit pass: `data/legacy` was restored from `data/legacy.bk` (byte-identical
copy), then every implemented `scripts/figures/*` script was executed
end-to-end against the live repo data. Where a script failed, the cause was
classified as either (a) an obvious, unambiguous authoring bug in the
already-verbatim-ported legacy source (typo, undefined variable, missing
`dir.create()`, disabled `ggsave()`) -- fixed here, minimally, with no
change to any computation, filter, threshold, column name, or file name --
or (b) missing input data, which is left unfixed per `CLAUDE.md`'s hard
constraint (report the gap, don't compute the data ourselves).

Every fix below is a single-token/single-line change (rename an undefined
variable to the one clearly intended by context, remove a stray trailing
comma, add a missing `dir.create()`, or uncomment an already-written
`ggsave()` call). No file name, column name, threshold, filter, or grouping
was invented or guessed anywhere -- where the fix would have required
assuming a file/column name that isn't already used elsewhere in the same
script, the script was left failing and is documented as blocked instead
(see "Blocked on missing data" below).

## Fixed: obvious authoring bugs (minimal, no logic change)

### `figureS1b_gleason_concordance.R`
- `ggsave(plot_path, p, ...)` referenced `save_dir` (never defined) and `p`
  (never assigned -- `(g_strip + g_main) + plot_layout(...)` was only
  auto-printed, not stored). Fixed: `save_dir` -> `figures_dir` (the only
  defined output-dir variable in the script); assigned the plot expression
  to `p` before `ggsave()`.

### `figureS3b_cluster_concordance.R`
- Same `save_dir` (undefined) -> `figures_dir` typo, one occurrence.

### `figureS3c_progression_km.R`
- `plot_path <- file.path(result_dir, plot_name)` referenced `result_dir`
  (never defined) -> `figures_dir`. The `ggsave(...)` call immediately below
  was also commented out (plot computed, never saved) -- uncommented it, no
  change to the plot object itself.

### `figureS5b_inflammation_km.R`
- `dir_inflam <- file.path(figures_dir, "inflammation")` was used as a
  `ggsave()` target but never `dir.create()`d (unlike its parent
  `figures_dir`, which is). Added the missing `dir.create(dir_inflam, ...)`
  call, mirroring the pattern already used for `figures_dir` two lines
  above.

### `figureS6a_stromogenic_km.R`
- Same missing-`dir.create()` pattern as `figureS5b_inflammation_km.R`, for
  `dir_stromo`.

### `figure6_inflammation_violin.R` (Figure 6d)
- `df_long$inflammation <- factor(df_long$inflammation, ...)` referenced
  `df_long`, which is never defined anywhere in the script -- the
  immediately-preceding dataframe (built by `df <- df_clr %>%
  inner_join(df_metadata, ...)`, which has the `inflammation` column) is
  unambiguously what was intended, since it's also exactly what
  `df_long_full` is derived from a few lines later via
  `df %>% pivot_longer(...)`. Fixed: `df_long` -> `df`, both occurrences.

### `figure7a_stromogenic_violin.R` (Figure 7a)
- `my_cols <- c(no = "#9efa70", yes = "#c55797", )` -- trailing comma inside
  `c(...)`, a hard R syntax error (unlike Python/JS, R does not allow
  trailing commas in argument lists). Removed the trailing comma.
- `ggplot(df_long_full, ...)` (x2) referenced `df_long_full`, which is never
  assigned -- `df_long` (built earlier via `df %>% pivot_longer(...)`, with
  the `niche`/`clr_proportion`/`direction` columns the plot needs) is what
  was intended. Fixed: `df_long_full` -> `df_long`, both occurrences.

### `figure7c_myCAF_km.R` (Figure 7c)
- `dir.create(figures, showWarnings = FALSE, recursive = TRUE)` referenced
  `figures` (never defined) instead of `figures_dir` (defined the line
  above). Fixed: `figures` -> `figures_dir`.

## Fixed 2026-08-05: yaml color-mapping type mismatch, on explicit instruction

### `figure6_niche_abundance_heatmap.R` (Fig 6a) -- yaml list vs. vector, FIXED
- Once its missing input (`niche_frequencies_per_tma_id.parquet`, resolved
  2026-08-05 -- see below) was present, the script failed:
  `Error: The color mapping should be a named vector or a function.`
  Root cause: `annotation_colors[[col]] <- custom_annotation_colors[[col]]`
  (line 90) assigned the raw `yaml::read_yaml()` result directly.
  `yaml::read_yaml()` parses every per-column color map in
  `resources/colormaps.yaml` (`pat_id`, `os_status`, etc.) as an R `list`,
  not a named character vector -- confirmed via
  `class(yaml::read_yaml("resources/colormaps.yaml")[["pat_id"]])` ->
  `"list"`. `ComplexHeatmap::rowAnnotation()`'s `col` argument requires a
  named vector or function, and rejects a list outright. Two lines below,
  the script's own `niche_colors <- unlist(colormap_niche)` shows the author
  was aware `yaml::read_yaml()` needs unlisting for this exact purpose, but
  never applied the same conversion to `annotation_colors[[col]]`.
  **Fixed per explicit user instruction** ("parse the yaml and convert into
  a format accepted by ComplexHeatmap, minimal changes only"): wrapped the
  assignment in `unlist()` --
  `annotation_colors[[col]] <- unlist(custom_annotation_colors[[col]])` --
  the exact same conversion the script already applies to `niche_colors` two
  lines below, just extended to this second call site. One line changed, no
  other logic touched.

### `figure6_niche_abundance_heatmap.R` (Fig 6a) -- `gleason_grp` label mismatch, FIXED 2026-08-05
- After the `unlist()` fix above, a second error surfaced:
  `Error: gleason_grp: cannot map colors to some of the levels: 1, 2, 5, 4, 3`.
  Root cause: `resources/colormaps.yaml`'s `gleason_grp` color keys are
  float-style strings (`"1.0"`, `"2.0"`, ..., `"5.0"`), while the script's
  own `hue_order_list` (line ~69) defines the factor levels as `"1"`,
  `"2"`, ..., `"5"` (no decimal) -- two independently-authored sources that
  disagree on label format. **Fixed per explicit user instruction**
  ("convert the yaml keys to ints to match the factor levels"): after the
  `unlist()`, for `gleason_grp` specifically, any numeric-looking key is
  reformatted via `as.character(as.integer(as.numeric(key)))` (`"1.0"` ->
  `"1"`); non-numeric keys (`"None"`, `"nan"`) are left untouched. Scoped to
  `gleason_grp` only -- no other column's keys are touched by this block.

### `figure6_niche_abundance_heatmap.R` (Fig 6a) -- `inflammation`/`stromogenic_smc_loss_reactive_stroma_present`/`glandular_atrophy_pin` yaml boolean-parsing gotcha, FIXED 2026-08-05
- Fixing `gleason_grp` surfaced a third, related error:
  `Error: inflammation: cannot map colors to some of the levels: no, yes`.
  Root cause: `resources/colormaps.yaml` writes these three columns' keys as
  bareword `no`/`yes` (unquoted), which YAML 1.1 (the spec `yaml::read_yaml()`
  implements) parses as **booleans**, not strings -- confirmed via
  `names(yaml::read_yaml("resources/colormaps.yaml")[["inflammation"]])` ->
  `"FALSE" "TRUE" "None" "nan"` instead of the intended `"no" "yes" "None"
  "nan"`. Same root cause, same class of bug as the `gleason_grp` fix above
  (a yaml-key-format vs. factor-level mismatch), so fixed the same way, per
  the same instruction: for any column whose parsed color-map names include
  both `"TRUE"` and `"FALSE"`, map `"FALSE"` -> `"no"` and `"TRUE"` ->
  `"yes"` (the only interpretation consistent with what
  `resources/colormaps.yaml` visibly intended to write). Applies generically
  in the same loop, not hardcoded to `"inflammation"` by name, so it also
  covers `stromogenic_smc_loss_reactive_stroma_present` (also in
  `df_metadata`) and `glandular_atrophy_pin` (not currently read by this
  script, but would hit the same issue if it were) without needing a
  separate fix each.
- With all three of these yaml/`ComplexHeatmap` mismatches fixed,
  `figure6_niche_abundance_heatmap.R` now runs end-to-end and saves
  `output/figures/figure6/figure6a_niche_proportion_heatmap_tma.pdf`.

## Resolved 2026-08-05: missing data supplied by Melissa

Melissa provided `data/melissa-transfer/{clinical,clusters_annotated_v2,
props_niche_tma_id}.parquet` and identified that the two previously-missing
assets are the same data under a different name, not separate files:
`cell_annotation.parquet` = `clusters_annotated_v2.parquet`,
`niche_frequencies_per_tma_id.parquet` = `props_niche_tma_id.parquet`.
Verified by `md5sum` -- both are **byte-identical** to the copies already
staged in `data/legacy/5-niches/`. Copied into place at the exact paths the
verbatim-ported scripts expect (filename only changes, no logic/column
changes): `data/cell_annotation.parquet` (from
`clusters_annotated_v2.parquet`) and
`data/legacy/5-niches/frequencies/niche_frequencies_per_tma_id.parquet`
(from `props_niche_tma_id.parquet`). `clinical.parquet` was also compared
cell-by-cell against `data/clinical.parquet` (after aligning on `sample_id`)
-- identical data, the only difference was the pyarrow/pandas writer
version (different file bytes, same content); no action needed.

This unblocked the five `cell_annotation.parquet`-dependent KM scripts,
which surfaced three further bugs, all now fixed (see below):

### `figure6_km_niche6.R`, `figure7b_niche9_km.R`, `figureS6bc_niche_km.R` -- missing `library(patchwork)`
- `p_prog$plot / p_prog$table` uses patchwork's `/` operator to combine two
  ggplot objects, but `library(patchwork)` was never loaded (unlike the
  sibling `figure7c_myCAF_km.R`, which already has it). Added the missing
  `library(patchwork)` call, matching the already-correct sibling.

### `figure6_km_niche6.R`, `figure7b_niche9_km.R`, `figureS6bc_niche_km.R`, `figure7c_myCAF_km.R` -- stale copy-pasted save gate
- The `ggsave()` for the combined plot+table panel was gated behind
  `if (col %in% c("luminal_infiltrated", "luminal_CAF1(CD105High)",
  "tumor_CAF1(CD105High)"))` (niche scripts) or
  `if (col == "epithelial-luminal(ERG+p53+)")` /
  `if (col == "stromal-CAF1(CD105+)")` (`figure7c_myCAF_km.R`) -- a
  hardcoded condition that was clearly copy-pasted from a shared/looped
  source and never updated after each script's `cols` was trimmed to its
  own single niche/label (per this project's established trimmed-loop
  porting pattern). Confirmed by execution: niche 6's panel and niche 8's
  panel (inside `figureS6bc_niche_km.R`) never saved a PDF at all before
  this fix, and `figure7c_myCAF_km.R`'s overall-survival panel (`p_os`)
  never saved either. Fixed by replacing the hardcoded condition with
  `col %in% cols` in all four scripts -- each panel now always saves for
  the niche/label it's actually built for, no other logic touched.
- Two of the three `ggsave()` calls this gate protects were also still
  commented out (`figure6_km_niche6.R`, `figure7b_niche9_km.R`,
  `figureS6bc_niche_km.R`) -- uncommented, same pattern as
  `figureS3c_progression_km.R` above.

### `figure6_km_niche6.R`, `figure7b_niche9_km.R`, `figureS6bc_niche_km.R`, `figure7c_myCAF_km.R` -- missing path separator
- `plot_name <- paste0(figures_dir, "km_survival_", ...)` concatenates
  without a `/`, e.g. producing
  `.../figure7/cell_typeskm_survival__disease_progr_...pdf` instead of
  `.../figure7/cell_types/km_survival__disease_progr_...pdf` -- confirmed
  by execution: `figure7c_myCAF_km.R`'s progression-free panel was actually
  written to `output/figures/figure7/` (one level up from its intended
  `cell_types/` subdirectory) under this malformed name. Fixed by wrapping
  the `paste0(...)` filename in `file.path(figures_dir, ...)` in all four
  affected scripts, matching the `file.path()` pattern already used
  elsewhere in the same scripts. The stray malformed file from the earlier
  run was deleted.

## Blocked on missing data (not fixed, per `CLAUDE.md`'s hard constraint)

None remaining as of 2026-08-05 -- both previously-missing assets
(`cell_annotation.parquet`, `niche_frequencies_per_tma_id.parquet`) were
supplied by Melissa and wired in, see above.

`figure6_niche_abundance_heatmap.R` (Fig 6a) is no longer blocked -- its
input data was supplied and all three yaml/`ComplexHeatmap` color-mapping
bugs are fixed, see above. It now runs end-to-end.

## Confirmed working, no bugs found

`figure2_cell_type_heatmap.R` (2a), `figure2_umap.py` (2b, S2a-c),
`figure3_caf_heatmap.R` (3b), `figure3_caf_umap.py` (3a),
`figure4_patient_clustering.py` (4a), `figure4_metagroup_barplot.py` (4b),
`figure4c_patient_cluster_km.R` (4c), `figure4de_cox_hazard_ratio.R` (4d-e),
`figure5_niche_heatmap.R` (5a), `figure5_niche_correlation.R` (5b),
`figure6c_niche_composition_filtered.py` (6c),
`figure7def_circos_plots.py` (7d-f, old `000_paper/11_niches` source --
see `open-questions.md`/`table.md` for the newer `sync_paper` source found
2026-08-05, not applied here without separate confirmation),
`figureS1a_cohort_summary.R` (S1a), `figureS3a_tma_stacked_barplot.py` (S3a),
`figureS4b_niche_mean_composition.py` (S4b),
`figure6_km_niche6.R` (6b), `figure6e_immune_risk_score_km.R` (6e),
`figure7b_niche9_km.R` (7b), `figure7c_myCAF_km.R` (7c),
`figureS6bc_niche_km.R` (S6b-c), `figure6_niche_abundance_heatmap.R` (6a) --
unblocked 2026-08-05 by Melissa's data and the yaml/`ComplexHeatmap` fixes,
see "Resolved 2026-08-05" above for the bugs found and fixed along the way.

## New 2026-08-05: `figure7_full_circos_plots.py` -- ported, blocked on a source bug

Per user request, `sync_paper/06-spatial-niches/interactions/interaction_compute_circos.py`
(the newer, full compute+plot circos source found earlier this session) was
ported verbatim to `scripts/figures/figure7_full_circos_plots.py`, kept
alongside (not replacing) the existing `figure7def_circos_plots.py` (older
`000_paper/11_niches` source, still works, see `figures.md`). Only
path/env-var substitutions and one disclosed environment-compatibility patch
were applied:
- `paths['clusters_annotated']` -> `LEGACY_DATA_DIR/5-niches/annotation/clusters_annotated.parquet`
- `paths['anndatas_dir']` -> the live legacy path (not staged, 476 pickles,
  too large to duplicate) --
  `/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas/`,
  confirmed matching legacy's own hardcoded path structure and 100%
  `sample_id` overlap with `clusters_annotated_v2.parquet` (see `data/assets.md`)
- `paths['colormaps']` -> `resources/colormaps.yaml`
- `paths['interaction_analysis_dir']` (data) / `paths_figures[...]` (figures)
  -> `$OUTPUT_FIGURES_DIR/figure7/full_circos/{data,}`
- **Disclosed environment-compatibility patch, not a logic change**: these
  anndata pickles were written by an older `anndata` whose
  `AnnDataFileManager` stored its backing-file reference under state key
  `_adata`; installed `anndata` 0.12.10 expects `_adata_ref` and raises
  `KeyError` otherwise (confirmed via direct unpickle test). Patched
  `AnnDataFileManager.__setstate__` to accept either key -- same category of
  fix as `scripts/port/port_umap_reducer.py`'s existing `pynndescent`
  compatibility patch. Confirmed these objects aren't actually file-backed
  (`_filename`/`_file` are `None`), so this is a pure key-rename shim.

**FIXED 2026-08-05, per explicit user instruction** ("if fixing this does
not change the logic please go ahead and retain the index"): originally
failed at `sample_list = df_clusters.index.get_level_values("sample_id")...`
with `KeyError: 'Requested level (sample_id) does not match index name
(None)'`. Root cause: the line immediately above,
`df_clusters.reset_index(inplace=True)`, converted
`clusters_annotated.parquet`'s `sample_id`/`object_id` MultiIndex into
regular columns -- but this line and the later `df_clusters.loc[sample]`
(per-sample lookup, further down) both require that MultiIndex to still be
in place. Fixed by removing the `reset_index(inplace=True)` call entirely,
retaining the MultiIndex `clusters_annotated.parquet` is already stored
with -- confirmed this is the correct, logic-preserving interpretation (not
a guess) because every downstream usage in the script (`.index.get_level_values`,
`.loc[sample]`, and `sample_data.index = sample_data.index.astype(str)`
followed by a positional merge with `adata.obs`) already assumes this exact
MultiIndex shape; no computation changed. Verified by execution: the script
now proceeds correctly past this point (`sample_data`/`obs` shapes match
for the first sample, `(3280, 8)` / `(3280, 10)`). Left running as a
long-running background job -- the per-edge interaction-counting loop
(`compute_interaction_matrices`) is unvectorized Python, iterating per
sample x per niche across all 476 samples, and is verbatim-slow by design
(not a bug, not optimized here to preserve the source's exact logic).

## OPEN ISSUE as of 2026-08-06: Fig 6a is currently BROKEN (fix reverted, not reapplied)

`figure6_niche_abundance_heatmap.R` was deliberately reverted to the pure
verbatim legacy version on 2026-08-06 (per explicit user instruction, to
diagnose the crash from first principles -- see below for exactly what
fails). It was never re-fixed afterward. **As currently committed, running
this script crashes** with `Error: The color mapping should be a named
vector or a function.` and produces no output.

The PDF at `output/figures/figure6/figure6a_niche_proportion_heatmap_tma.pdf`
is **stale** -- it's the output of the fixed version described below,
generated before the revert. It still visually matches the published
panel, but the current script state cannot reproduce it. The fix (keep
`tma_id` for indexing only, drop it from the displayed row annotation, see
below) needs to be re-applied.

## New 2026-08-05: Fig 6a `tma_id` fix was incomplete -- fixed properly (NOT CURRENTLY APPLIED, see OPEN ISSUE above)

Earlier this session, `tma_id` was uncommented back into `select()` to fix
a genuine runtime crash (`df_metadata$tma_id` didn't exist, producing a
0-row matrix and a downstream `ComplexHeatmap` color-mapping error). That
fix was directionally necessary but incomplete: it silently changed TWO
things, not one. (1) It correctly restored the full 459-row TMA granularity
-- confirmed empirically: `distinct()` on the clinical columns *without*
`tma_id` collapses 459 rows down to 346, because many patients have
multiple TMA cores sharing identical clinical field values. (2) As an
unexamined side effect, `tma_id` also became a visible row-annotation
track, since `rowAnnotation(df = df_metadata, ...)` displays every column
passed to it -- and direct visual comparison against the published Figure
6a (rendered via `pymupdf`, installed this session for this purpose)
confirmed the paper shows no `tma_id` annotation track, only
pat_id/os/disease_progr/gleason_grp/inflammation/stromogenic.

Fixed properly: `tma_id` stays in `df_metadata` (needed for `distinct()`
and row-reordering), but a separate `df_display <- df_metadata %>%
select(-tma_id)` is passed to `rowAnnotation()` instead -- keeps the
correct 459-row data, drops only the display column. `show_legend` updated
accordingly (`c(FALSE, rep(TRUE, ncol(df_display) - 1))` -- hide only
pat_id's legend, matching the paper's legend block of exactly 5 categorical
variables). Re-rendered and visually compared against the paper: row count
and overall dendrogram/cluster block structure now closely match (same 5
main clusters, same general niche groupings) -- the earlier "clustering
doesn't look the same as the paper's" discrepancy was very likely driven
by this 459-vs-346-row skew, not a distance-function or algorithm
difference (both were already byte-identical to the legacy source).

## New 2026-08-05: S2a-c was wrongly documented as already covered

`figures.md` previously claimed Supplementary Fig 2a-c was already produced
by `figure2_umap.py`'s per-`main_group` subset loop (colors/filters points
from the single all-cells UMAP embedding). While checking for other
compartment-specific `reducer.pkl` files (per user request, immune/
epithelial/endothelial), found a third legacy script,
`archive/scripts/02-umaps/0-umaps-main-types.py`, that fits a genuinely
**separate UMAP per compartment** -- its own `reducer.pkl` per
`main_group` at `/work/.../data/PCa/0-paper/2-umaps/2-main_groups/`
(confirmed 5 exist: immune 677MB, epithelial 2.2GB, endothelial 168MB,
stromal 1.4GB, undefined 122MB -- only immune/epithelial/endothelial are
part of S2 per its legend). This means the previous S2a-c output did NOT
match the published supplementary figure (wrong embedding entirely, not
just a cosmetic difference).

Per explicit user instruction, ported all three relevant `reducer.pkl`
files (same `port_umap_reducer.py` mechanism, same pinned pixi env) to
`data/figures/figureS2_main_groups_umap/{immune,epithelial,endothelial}/umap_embeddings.parquet`,
and wrote a new dedicated `scripts/figures/figureS2_compartment_umap.py`
(1:1 port of `0-umaps-main-types.py`'s `label`/per-marker plotting loops,
does not touch `figure2_umap.py`). `figures.md`'s S2 row corrected.

## New 2026-08-05: Fig2b/3a UMAP embeddings re-ported from scratch, per explicit user instruction

Earlier this session, Fig2b/3a's embeddings were restored from
`output/figures.bk2/` (a prior session's already-computed copy) after an
unwanted from-scratch UMAP *re-fit* was killed mid-run (see below). Per
explicit user instruction ("I don't want to reuse the previously ported
reducer.pkl but port them from scratch"), re-ran the actual port from the
raw legacy `reducer.pkl` files via `port_umap_reducer.py`, confirming the
correct source paths against `archive/scripts/02-umaps/0-umaps.py` (Fig 2)
and `0-umaps-cafs.py` (Fig 3)'s own `params`/`get_reducer_path()` logic --
not guessed. Output: `data/figures/figure2_umap/umap_embeddings.parquet`,
`data/figures/figure3_caf_umap/{excl_markers,caf_markers_only}/umap_embeddings.parquet`.
Values confirmed identical to the `.bk2` copy (same underlying `reducer.pkl`,
different derivation path), but this is now a clean, from-scratch,
independently-verified re-derivation rather than a reused prior artifact.
Documented the exact `reducer.pkl` paths and rerun commands directly in
`scripts/data/figure2_umap_embedding.py` and
`scripts/data/figure3_caf_umap_embedding.py`'s docstrings, per explicit
user instruction, so this port can be redone without re-deriving the paths
from the archive scripts each time.

**Important correction, also per explicit user instruction**: earlier in
this session, when `figure2_umap.py`/`figure3_caf_umap.py` failed because
`output/figures/` (containing their cached embeddings) had been moved
aside, the response was to kick off `scripts/data/figure2_umap_embedding.py`/
`figure3_caf_umap_embedding.py` -- which do NOT port anything, they fit a
**fresh, unseeded UMAP** each time (their own docstrings, added this
session, now state this explicitly: "NOT reproducible against the
published figure"). This was wrong and was killed before completion, per
direct user correction ("We want to port them with the porting scripts not
recompute from scratch"). The two `_embedding.py` scripts are left as-is
(unchanged logic) but now carry a prominent docstring warning against
using them as the source of truth for the published figures.

## New 2026-08-05: `figureS7_ari_robustness.py` ported and REPRODUCED

Per user request, ported `sync_paper/06-spatial-niches/construction/kmeans_clustering.py`'s
ARI-robustness section (lines ~217-282 of that file) to
`scripts/figures/figureS7_ari_robustness.py`. Reuses the exact same data
loading and `perform_kmeans_clustering`/`wrapper_nhood_filtering` utilities
as the already-validated `figure5_niche_clustering.py` (same
`LEGACY_DATA_DIR/PCa_NHood/CellCellNeighborhoods/` source, same
`PCA_NHOODs_clean/robustness/utils/` imports) -- does not modify
`figure5_niche_clustering.py`. Verbatim port: `k=24`, `n_runs=50`, seeds via
`np.random.RandomState(42)`, `construct_ari_matrix()`/`prepare_ari_data()`,
boxplot of pairwise ARI per run with the best-agreement run (mean ARI)
highlighted yellow -- confirmed this is `seed=686`, the same seed
`figure5_niche_clustering.py` already uses for the real clustering,
matching the legacy source's own comment ("seed with top ARI mean across
runs"). No bugs found in this section during the port.

Ran as a long-running background job: 2,051,915-row neighborhood graph x 50
separate k-means(k=24) fits, each followed by ARI computation against all
previous runs (O(n_runs^2) pairwise ARI). **Completed successfully**
(~2h15m runtime), no errors -- only a harmless seaborn deprecation warning
(`palette` without `hue`, cosmetic, not a bug). Output saved to
`output/figures/figureS7/ari_boxplot.png`. S7 is confirmed reproducible.

## New 2026-08-05: Fig 7b (niche 9) REPRODUCED -- visually confirmed against the paper

`figure7b_niche9_km.R`'s output
(`figure7/niches/km_survival__disease_progr_luminal_CAF1(CD105High)with_table.pdf`)
was directly compared against the published Fig 7b panel (paper page
rendered via `pymupdf`). Number-at-risk table matches almost exactly:
paper shows low 75/54/20/0, high 115/58/13/0 at t=0/50/100/150; our output
shows the identical 75/54/20/0 and 115/58/13/0. Displayed p-value
(raw 0.0019 vs. paper's 0.018) has the same explanation as Fig 6b/S6b/S6c/
7c below -- the legacy source always plots the raw log-rank p-value, never
the BH-adjusted `qval` it separately computes across the full niche family.
Confirmed reproduced.

## New 2026-08-05: Fig 7c (myCAF) p-value mismatch -- same adjustment issue as Fig 7b/S6b/S6c

`figure7c_myCAF_km.R`'s displayed progression-free p-value (0.0137, matches
legacy `celltype_km.R`'s always-raw `ggsurvplot(pval = TRUE, ...)`) doesn't
match the paper's reported p = 0.24 for the same panel. Root cause
identical to the Fig 7b/S6b/S6c finding above: legacy `celltype_km.R`
(un-trimmed) loops over all 35 labels (`cols <-
colnames(df_props)[2:ncol(df_props)]`) and computes a BH-adjusted `qval`
across that full family, written only to a CSV, never back into the plot.
Verified empirically: re-running the full 35-label family and BH-adjusting
jointly gives `stromal-CAF1(CD105+)`'s qval = 0.2389, matching the paper's
0.24 almost exactly. Not a bug in our port; same accepted-as-is status as
the niche KM panels. Also confirmed the myCAF label identity is correct
(not a guess): the paper's own text (page 8) states "CD105high are
annotated as myCAFs," matching `cols <- "stromal-CAF1(CD105+)"` exactly.

## New 2026-08-05: S4b orientation differs from paper, not a data bug

Spot-checked several niches' segment proportions in
`figureS4b_niche_mean_composition.py`'s output against the published panel
(paper page rendered via `pymupdf`) -- they match closely. The visible
difference is orientation: the paper's panel is a horizontal stacked bar
(niches as rows), ours is vertical (niches as columns), inherited directly
from the legacy source
(`000_paper/11_niches/111_heatmaps/visualize_composition.py`'s
`df_comp_mean.plot(kind='bar', ...)` -- also vertical). Our script is a
faithful verbatim port; the paper's horizontal orientation isn't produced
by any script found in the legacy repo -- likely a manual
post-processing/rotation step before publication, same category as the
p-value substitution found for the KM panels above. Not fixed (would be a
disclosed cosmetic deviation from the verbatim `kind='bar'` call).

## Not run (no implemented script)

Fig 1a-c, 3c-e, 5c, 7g-h, S4a -- non-code-derived or no legacy plotting
code exists (see `figures.md`).

`figure5_niche_clustering.py` and `figure5_niche_annotation.py` were not
re-run in this pass -- they're data-preparation scripts (not plotting
scripts), and their output (`clusters_annotated_v2.parquet`,
`niche_annotations_v2.csv`) is already staged and verified byte-identical
in `data/legacy/5-niches/annotation/` (see `data/assets.md`).
