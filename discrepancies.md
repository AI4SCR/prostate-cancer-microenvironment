# Logic-fidelity audit: `scripts/figures/*` vs. legacy (old repo) sources

Each `scripts/figures/figureN_*` file is supposed to be a faithful port of an
old-repo script (`000_paper/...`), changing **only** hardcoded paths,
uninstallable-package substitutions (documented), and cosmetic cleanup — no
computational logic. This audit (11 parallel sub-agents, one per script pair)
checked that claim file by file. Status below reflects findings at audit
time; ✅ = fixed since, ⬜ = not yet fixed.

## Summary

| Script | Verdict | Real discrepancies | Status |
|---|---|---|---|
| `figure5_niche_correlation.R` | PASS | 0 | — |
| `figure6_km_niche6.R` | PASS | 0 (disclosed CSV additions only) | — |
| `figure4_survival.R` | PASS | 0 | — |
| `figure5_niche_clustering.py` | FAIL | 2 minor | ⬜ |
| `figure5_niche_annotation.py` | FAIL | 2 minor + 1 trivial | ⬜ |
| `figure5_niche_heatmap.R` | FAIL | 2 real | ⬜ |
| `figure6_niche_abundance_heatmap.R` | FAIL | 2 real, 1 cosmetic | ⬜ |
| `figure6_inflammation_violin.R` | FAIL | 1 real (missing layer, 2 plots) | ⬜ |
| `figure2_umap.py` | **FAIL** | **6 major** | ⬜ |
| `figure2_cell_type_heatmap.R` | **FAIL** | **8 major** | ⬜ |
| `figure3_caf_umap.py` | **FAIL** | **3 major** | ⬜ |

The Figure 2/3 failures are the most serious: those scripts were reported
as "validated against real data" in earlier work on this repo, but
"validated" only meant *ran without error and produced a plot* — nobody had
checked the actual computation against the old repo's script line-by-line
until now. They do not currently reproduce the paper's method.

---

## `figure2_umap.py` — 6 major discrepancies

vs. `000_paper/02_umaps/0-umaps.py` (+ `0-umaps-main-types.py`, not otherwise used)

1. **UMAP `n_neighbors`: 50 → 15.** Legacy passes `50`; new script defaults to `15`. Changes the embedding.
2. **Normalization step removed entirely.** Legacy calls `normalize(data, exclude_zeros=True)` before UMAP; new script has no equivalent call.
3. **Marker exclusion removed.** Legacy excludes `{dna1, dna2, icsk1, icsk2, icsk3, fap}` before computing the embedding; new script uses `marker_columns(cells)` with no such exclusion.
4. **`metric="euclidean"` added.** Not present in the legacy call; UMAP's own default differs from an explicit euclidean setting in some configurations.
5. **Marker plotting scope narrowed.** Legacy plots every marker column; new script hardcodes 4 markers (`cd45`, `pan_keratin`, `cd31`, `vimentin`).
6. **Data source changed.** Legacy loads via `PCa()` + `normalize()`; new script uses `load_exported_cells(export_dir, exclude_undefined=False)` — plausibly equivalent in principle but not verified equivalent, and moot given (2) above (the normalization this alternate path would need is simply missing).

## `figure2_cell_type_heatmap.R` — 8 major discrepancies

vs. `000_paper/04_heatmaps/2-cell-types-heatmap.R` (the `heatmap()`/`heatmap.agg()` functions)

1. **Filter target differs.** New filters `label != "undefined"` (a disclosed, git-logged fix); legacy's `heatmap.agg()` instead filters out `label == 'mix-vessels-PMN-MDSCs'` — these are different cell populations, not the same fix expressed differently.
2. **No subsampling.** Legacy subsamples to 2000 cells when not aggregating; new script has no sampling step (though it does always aggregate, see #7, which may make this moot — needs a decision, not an assumption).
3. **Matrix orientation flipped.** Legacy transposes so markers are rows; new script keeps cell types as rows, markers as columns.
4. **Explicit clustering method/distance added.** New hardcodes `clustering_method_rows = "average"`, `clustering_distance_rows = "euclidean"`; legacy leaves these at ComplexHeatmap defaults (which may differ).
5. **Fewer annotations.** New shows only `main_group` as a `right_annotation`; legacy's `heatmap()` shows `label`, `main_group`, AND `patient` as `top_annotation`.
6. **Color source changed.** New generates `main_group` colors via `circlize::rand_color()` (non-deterministic across runs unless seeded); legacy loads fixed colors from `colormaps.yaml`.
7. **Aggregation always on.** New always does `group_by(label)`; legacy's `heatmap()` only aggregates conditionally (`aggregate_by != NULL`), otherwise plots individual subsampled cells.
8. **New output added.** New writes `figure2a_mean_expression.parquet`, which legacy never did — likely fine as a disclosed addition, but wasn't disclosed as such.

## `figure3_caf_umap.py` — 3 major discrepancies

vs. `000_paper/02_umaps/0-umaps-cafs.py`

1. **Cell filter differs.** New: `cells["main_group"] == "stromal"` (all stromal cells). Legacy: `metadata.label.str.contains('CAF')` (only labels containing "CAF" — a narrower, different set; stromal includes non-CAF stromal cells like pericytes).
2. **UMAP `n_neighbors`: 50 → 15.** Same issue as `figure2_umap.py`.
3. **Normalization removed.** Legacy calls `normalize(df, exclude_zeros=True)` before `compute_umap()`; new script feeds raw marker values directly to UMAP.

---

## Minor/real issues in Figure 5/6 scripts

### `figure5_niche_clustering.py`
- Output PNG filename hardcoded (`figure5_niche_kmeans_raw_heatmap.png`) instead of derived from `final_cluster_name` like legacy's `f"{nhood}_heatmap.png"`. Cosmetic but not disclosed.
- `engine="fastparquet"` dropped from the final `to_parquet()` call (legacy specifies it explicitly).

### `figure5_niche_annotation.py`
- `engine="fastparquet"` dropped from both `read_parquet()` and `to_parquet()` calls.
- Legacy repeats the exact same niche-mapping assignment twice (lines 56-58 and again 134-136, byte-identical, clearly redundant); new script does it once. Almost certainly harmless (idempotent), but is technically a removed step under a strict reading.

### `figure5_niche_heatmap.R`
- `my_colors`/`col_fun_presence` computation adds `na.rm = TRUE` to `range(presence, ...)` — legacy's `range(presence)` has no NA handling. Only matters if `presence` actually contains NAs.
- New script's final `draw()` call includes `annotation_legend_list = list(lgd_presence)`; legacy's real (uncommented, inside-pdf) `draw()` call omits it — legacy has two earlier debug `draw()` calls (before opening the pdf device) that DO include it, so the omission in legacy's actual saved output looks like a legacy bug, not intent. New script's behavior matches legacy's evident intent, not its literal final line.

### `figure6_niche_abundance_heatmap.R`
- Legacy's heatmap `pdf()`/`dev.off()` calls are commented out — the legacy script, as written, never actually saves this heatmap to disk. New script enables saving. Given the entire point of a figure-reproduction script is to produce the figure, this is very likely a legacy authoring artifact (script left mid-edit) rather than intentional suppression — but it's still a literal behavior change from "no file written" to "file written."
- Progression-plot y-axis label: legacy says `"Survival probability"` for the progression panel too (an apparent legacy copy-paste bug, since it's plotting `clinical_progr`); new script corrects it to `"Progression-free survival probability"`.
- Output filename changed (organizational, matches this repo's `figureN*` naming convention elsewhere — low concern).

### `figure6_inflammation_violin.R`
- **Real, unexplained omission**: legacy's `p_up` and `p_split` plots both include a `stat_summary(fun.data = "mean_se", geom = "pointrange", ...)` layer (mean ± SE marker on top of the violin/boxplot); the new script's equivalent two plots drop this layer entirely, even though the new script's own header comment claims "same statistics." The third plot (`p_sep`) does keep it. This is a straightforward bug in the port, not a considered substitution.

---

## Fix plan

All of the above will be corrected to match the legacy scripts exactly,
except where a deviation is already required and disclosed (survminer →
ggsurvfit, introdataviz → geom_violin+dodge, missing packages). Tracked and
applied in the commits immediately following this document.
