# Figure → script mapping

# Overview

Every panel of every figure in the paper, the script that produces it, its
legacy source, the input data it reads, and its output file(s). `Validated`
= directly visually confirmed against the published panel (`paper.pdf`);
`false` covers both known mismatches (see the Issue sections below) and
panels never explicitly checked against the paper. `—` = not applicable
(no script exists, either non-code-derived or a confirmed gap, see Issues).
All `sync_paper` paths are relative to `000_paper/sync_paper/`; other legacy
paths relative to `PCa/`.

| Figure | Panel | Script | Ported from | Assets | Output | Validated |
|---|---|---|---|---|---|---|
| 1 | a-c | — | — (BioRender + raw image crops) | — | — | — |
| 2 | a | `figure2_cell_type_heatmap.R` | `000_paper/04_heatmaps/2-cell-types-heatmap.R` | intensity_normalized.parquet, metadata.parquet | `figure2/figure2a_cell_type_heatmap.pdf` | false |
| 2 | b | `figure2_umap.py` | `000_paper/02_umaps/0-umaps.py` | metadata.parquet, intensity_normalized.parquet, clinical.parquet, `data/figures/figure2_umap/umap_embeddings.parquet` (ported from legacy `reducer.pkl`) | `figure2/label=*.pdf`, `figure2/value=*.pdf` (44 files) | true |
| 3 | a | `figure3_caf_umap.py` | `000_paper/02_umaps/0-umaps-cafs.py` | metadata.parquet, intensity.parquet, clinical.parquet, `data/figures/figure3_caf_umap/{excl_markers,caf_markers_only}/umap_embeddings.parquet` (ported from legacy `reducer.pkl`) | `figure3/{excl_markers,caf_markers_only}/label=*.pdf`, `value=*.pdf`, `cell_type=*.pdf` | false |
| 3 | b | `figure3_caf_heatmap.R` | `000_paper/04_heatmaps/2-cell-types-heatmap.R` (`heatmap.caf()`) | intensity_normalized.parquet, metadata.parquet | `figure3/figure3b_caf_heatmap.pdf` | false |
| 3 | c-e | — | — (representative ROIs, not code-derived) | — | — | — |
| 4 | a | `figure4_patient_clustering.py` | `sync_paper/05-heterogeneity/stacked-frequencies-label.py` (`pat_id` branch) | metadata.parquet, clinical.parquet | `figure4/figure4a_stacked_barplot.pdf`, `figure4a_patient_composition.parquet` | true |
| 4 | b | `figure4_metagroup_barplot.py` | `sync_paper/05-heterogeneity/stacked-frequencies-label.py` (tail block) | `figure4/figure4a_patient_composition.parquet` (Fig 4a's own output), resources/metalabels.yaml | `figure4/figure4b_metagroup_barplot.pdf` | true |
| 4 | c | `figure4c_patient_cluster_km.R` | `sync_paper/03_survival/patient_risk_group_km.R` | `5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet`, clinical.parquet | `figure4/figure4c_survival_by_patient_cluster.png` | true |
| 4 | d-e | `figure4de_cox_hazard_ratio.R` | `000_paper/03_survival/survival-proportions.r` | metadata.parquet, clinical.parquet | `figure4/figure4_cox_d.png`, `figure4_cox_e.png` | true |
| 5 | a | `figure5_niche_heatmap.R` (+ `figure5_niche_clustering.py`, `figure5_niche_annotation.py`) | `sync_paper/06-spatial-niches/composition/z_score_composition_vis.R`; `sync_paper/.../construction/kmeans_clustering.py`; `000_paper/11_niches/110_analysis/01_annotation_v2.py` | `5-niches/annotation/clusters_annotated_v2.parquet`, `niche_annotations_v2.csv`, `visualization/composition/niche_heatmap_data.parquet` | `figure5/figure5a_niche_zscore_heatmap.pdf` | false |
| 5 | b | `figure5_niche_correlation.R` | `sync_paper/06-spatial-niches/abundance/correlation_frequencies.R` | `5-niches/frequencies/stacked_barplots/props_niche_tma_id.parquet`, `annotation/niche_annotations_v2.csv`, clinical.parquet | `figure5/figure5b_niche_correlation_heatmap.pdf` | false |
| 5 | c | — | — (no legacy plotting code for per-core composition bars) | — | — | — |
| 6 | a | `figure6_niche_abundance_heatmap.R` | `sync_paper/06-spatial-niches/abundance/heatmap_frequencies.R` | `5-niches/frequencies/niche_frequencies_per_tma_id.parquet`, clinical.parquet | `figure6/figure6a_niche_proportion_heatmap_tma.pdf` (**stale**, see Issue) | false |
| 6 | b | `figure6_km_niche6.R` | `sync_paper/03_survival/niche_km.R` | clinical.parquet, cell_annotation.parquet | `figure6/niches/km_survival__disease_progr_tumorERG+p53+_ProlifLuminalwith_table.pdf` | true |
| 6 | c | `figure6c_niche_composition_filtered.py` | `sync_paper/06-spatial-niches/composition/niche_composition.py` (niches 6/16/17/18) | `5-niches/annotation/clusters_annotated_v2.parquet` | `figure6/figure6c_niche_composition_barplot.pdf` | false |
| 6 | d | `figure6_inflammation_violin.R` | `sync_paper/06-spatial-niches/histology/inflammation_vis.R` | `5-niches/annotation/clusters_annotated_v2.parquet`, clinical.parquet | `figure6/violin_boxplot_inflammation.pdf` | true |
| 6 | e | `figure6e_immune_risk_score_km.R` | `sync_paper/03_survival/inflammatory_niche_risk_group.R` | clinical.parquet, cell_annotation.parquet | `figure6/kaplan_meier_inflammation_{os,progression}_risk_group_max.pdf` | false |
| 7 | a | `figure7a_stromogenic_violin.R` | `sync_paper/06-spatial-niches/histology/stromogenic_vis.R` | `5-niches/annotation/clusters_annotated_v2.parquet`, clinical.parquet | `figure7/violin_boxplot_stromogenic.pdf`, `violin_boxplot_stromogenic_sep.pdf` | false |
| 7 | b | `figure7b_niche9_km.R` | `sync_paper/03_survival/niche_km.R` (niche 9) | clinical.parquet, cell_annotation.parquet | `figure7/niches/km_survival__disease_progr_luminal_CAF1(CD105High)with_table.pdf` | true |
| 7 | c | `figure7c_myCAF_km.R` | `sync_paper/03_survival/celltype_km.R` (`stromal-CAF1(CD105+)` = myCAF) | clinical.parquet, cell_annotation.parquet | `figure7/cell_types/km_survival__disease_progr_stromal-CAF1(CD105+)with_table.pdf` | true |
| 7 | d-f | `figure7def_circos_plots.py` | `000_paper/11_niches/114_interactions/visualize_circos_plot.py` (+ `circos_plots.py`) | `5-niches/visualization/interactions/.../dataframes_v2/{niche}.parquet` | `figure7/figure7d_luminal_infiltrated_circos_plot.pdf`, `figure7e_tumor_CAF1_lymphocytes_circos_plot.pdf`, `figure7f_luminal_CAF1(CD105High)_circos_plot.pdf` | true |
| 7 | g-h | — | — (representative ROIs, not code-derived) | — | — | — |
| S1 | a | `figureS1a_cohort_summary.R` | `000_paper/01_clinical_metadata/clinical.r` | clinical.parquet | `figureS1/{tma-level,patient-level}/*.pdf` (36 files) | true |
| S1 | b | `figureS1b_gleason_concordance.R` | `sync_paper/05-heterogeneity/patient-core-heterogeneity.R` (Gleason section) | clinical.parquet | `figureS1/heatmap_gleason_grp_by_gs_grp.pdf` | true |
| S2 | a-c | `figureS2_compartment_umap.py` | `000_paper/02_umaps/0-umaps-main-types.py` | metadata.parquet, intensity_normalized.parquet, `data/figures/figureS2_main_groups_umap/{immune,epithelial,endothelial}/umap_embeddings.parquet` (ported from legacy per-compartment `reducer.pkl`) | `figureS2/{immune,epithelial,endothelial}/label=label.pdf`, `value=*.pdf` | true |
| S3 | a | `figureS3a_tma_stacked_barplot.py` | `sync_paper/05-heterogeneity/stacked-frequencies-label.py` (`tma_id` branch) | metadata.parquet, clinical.parquet | `figureS3/figureS3a_stacked_barplot.pdf` | true |
| S3 | b | `figureS3b_cluster_concordance.R` | `sync_paper/05-heterogeneity/patient-core-heterogeneity.R` (cluster section) | `5-niches/barplot_data/metadata_with_dendrogram_colors_label_{pat_id,tma_id}.parquet` | `figureS3/heatmap_cluster_group_tma_by_patient_cluster_group.pdf` | true |
| S3 | c | `figureS3c_progression_km.R` | `sync_paper/03_survival/patient_risk_group_km.R` | `5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet`, clinical.parquet | `figureS3/progr_patient_cluster_group_final_all.pdf` | true |
| S4 | a | — | — (representative niche images, not code-derived) | — | — | — |
| S4 | b | `figureS4b_niche_mean_composition.py` | `sync_paper/06-spatial-niches/composition/niche_composition.py` | `5-niches/annotation/clusters_annotated_v2.parquet` | `figureS4/figureS4b_niche_composition_barplot.pdf` | false |
| S5 | a | — | — (no legacy plotting code for per-core composition bars) | — | — | — |
| S5 | b | `figureS5b_inflammation_km.R` | `sync_paper/03_survival/clinical_risk_group_km.R` (inflammation section) | clinical.parquet | `figureS5/inflammation/progr_inflammation.pdf`, `survival_inflammation_with_table.pdf` | true |
| S6 | a | `figureS6a_stromogenic_km.R` | `sync_paper/03_survival/clinical_risk_group_km.R` (stromogenic section) | clinical.parquet | `figureS6/stromogenic/progr_stromogenic_with_table.pdf`, `survival_stromogenic_with_table.pdf` | true |
| S6 | b-c | `figureS6bc_niche_km.R` | `sync_paper/03_survival/niche_km.R` (niches 2/8) | clinical.parquet, cell_annotation.parquet | `figureS6/niches/km_survival__disease_progr_luminal_infiltratedwith_table.pdf` (S6b), `..._tumor_CAF1_lymphocyteswith_table.pdf` (S6c) | true |
| S7 | a | `figureS7_ari_robustness.py` | `sync_paper/06-spatial-niches/construction/kmeans_clustering.py` (ARI robustness section) | `PCa_NHood/CellCellNeighborhoods/{cell_metadata,metadata}.parquet`, `graph_type=radius-radius=32/data.parquet` | `figureS7/ari_boxplot.png` | true |

Assets are relative to `LEGACY_DATA_DIR` unless otherwise noted; `metadata.parquet`/`clinical.parquet`/`intensity(_normalized).parquet`/`cell_annotation.parquet` are relative to `EXPORT_DIR`. A second script, `figure7_full_circos_plots.py` (ported from the newer `sync_paper/06-spatial-niches/interactions/interaction_compute_circos.py`), exists alongside `figure7def_circos_plots.py` for Fig 7d-f — see Issue below.

# Issue: Figure 6a currently broken

`figure6_niche_abundance_heatmap.R` was reverted to the pure verbatim
legacy version (to diagnose a crash from first principles) and was never
re-fixed. As currently committed, `tma_id` is commented out of `select()`
and the script crashes with `Error: The color mapping should be a named
vector or a function.` (root cause: the missing `tma_id` produces a 0-row
matrix; a separate `yaml`-list-vs-vector bug in the color mapping hits
first). The PDF on disk is stale, left over from a working, fixed version.

The verified fix (not currently applied): keep `tma_id` in `df_metadata`
for `distinct()`/row-indexing (dropping it collapses the correct 459 TMA
rows down to 346, since several patients share identical clinical values
across multiple cores), but build a separate `tma_id`-free data frame for
`rowAnnotation()` — the paper's published panel shows no `tma_id`
annotation track. With this fix, row count and dendrogram/cluster
structure matched the paper closely (same 5 main clusters).

# Issue: Figure 3b CAF heatmap doesn't match the paper

Pericytes aren't colored yellow, rows aren't clustered, and the overall
heatmap looks structurally different from the published panel — the
biggest visual mismatch found. Contributing factors already identified but
not fully resolved:
- No legacy script produces the paper's exact 12-marker panel from code
  alone (`2-cell-types-heatmap.R`'s `heatmap.caf()`: 9 markers;
  `2-1-cell-types-heatmap.R`'s newer sibling: 13 markers with `pdpn`
  commented out). The current script hardcodes a 12-marker list inferred
  by taking the 13-marker list, uncommenting `pdpn`, and dropping
  `c_casp3`/`ki_67` — reconstructed by inference + visual confirmation of
  the marker set only, not derived from a script.
- The two legacy `heatmap.caf()` versions also disagree on color scale
  (`lightblue/white/lightcoral` vs `#2166ac/white/#b2182b`); current script
  uses the older scale, unconfirmed against the published figure.
- Row clustering and pericyte color are not yet explained.

# Issue: Figure 7a — unclear which script actually produces this panel

Multiple scripts in the legacy repo produce violin plots in this area of
the codebase (stromogenic status vs. niche abundance); it isn't confirmed
which one is the actual source of the published Fig 7a panel. Current
script ports `sync_paper/06-spatial-niches/histology/stromogenic_vis.R`
per the general "always cite the newest `sync_paper` sibling" policy, but
this hasn't been independently verified against the published figure the
way other panels have.

# Issue: Supplementary Figure 4b orientation differs from the paper

Segment proportions match the paper closely when spot-checked niche by
niche, but the paper's panel is a horizontal stacked bar (niches as rows)
while ours (and the legacy source itself — confirmed via
`000_paper/11_niches/111_heatmaps/visualize_composition.py`'s own
`df_comp_mean.plot(kind='bar', ...)`) is vertical. Our script is a
faithful verbatim port; the paper's horizontal orientation isn't produced
by any script found — likely a manual rotation before publication, never
captured in code. Not fixed (would be a disclosed cosmetic deviation from
the verbatim `kind='bar'` call).

# Issue: Figure 5a/5b — verbatim port doesn't match the paper's layout

Both scripts are confirmed byte-for-byte verbatim ports (checked against
all legacy copies of each source), so these are documented, unexplained
mismatches, not porting bugs:
- **5a**: legacy code puts the column annotation (`top_annotation`) at the
  top of the heatmap in every copy found; the published figure has it at
  the bottom.
- **5b**: legacy code draws per-cell correlation numbers via `cell_fun`
  and lets `ComplexHeatmap` choose row/column order; the published figure
  has no numbers in the cells and a different (mirrored) row/column order.

# Issue: Figure 6a/6b/7b/7c/S6b/S6c — paper displays an adjusted p-value the plotted script never computes

Confirmed mechanism, not a bug: the legacy KM scripts (`niche_km.R`,
`celltype_km.R`, un-trimmed) always plot the raw, unadjusted log-rank
p-value via `ggsurvplot(pval = TRUE)`. Separately, after looping over
every niche/label, they compute a BH-adjusted `qval` across that full
family — written only to a CSV, never fed back into any plot. Verified
empirically for three cases: re-running the full un-trimmed family and
BH-adjusting jointly reproduces the paper's displayed p-value almost
exactly (niche 2 progression-free: raw 0.393 → adjusted 0.6812 ≈ paper's
0.68; niche 8: same adjusted value, 0.6812; myCAF progression-free: raw
0.0137 → adjusted 0.2389 ≈ paper's 0.24). The paper's authors evidently
substituted the adjusted value into the published panel by hand. Not a bug
in any of these scripts; accepted as-is — the survival curves and
number-at-risk tables themselves match the paper (confirmed for 6b/7b/7c).

# Issue: Figure 5c / S5a — no legacy code for per-core composition bars

Both legends describe a per-representative-core, x=niche composition bar
(5c: "stacked barplots show the mean cell type composition of the
dominant niches in these cores"; S5a: "ordered by its proportion in each
core"). Re-searched `sync_paper/06-spatial-niches/`, `000_paper/11_niches/`,
and the older archive tree (including
`111_heatmaps/visualize_composition.py`, not previously checked) — no
script anywhere plots per-core composition this way, only the corpus-wide
mean-per-niche composition already ported as `figureS4b`/`figure6c`.
Figure 6c's own "two example cores" bottom panel has the same gap. Not
attempted, per this project's rule against writing new, non-ported
plotting code.

# Issue: Figure 7d-f — two circos scripts, the newer one is unfinished

`figure7def_circos_plots.py` (working, plots precomputed
`dataframes_v2/*.parquet`) ports the older `000_paper/11_niches` source.
A newer `sync_paper/06-spatial-niches/interactions/interaction_compute_circos.py`
was found and ported separately to `figure7_full_circos_plots.py` — it
computes interaction matrices itself from raw anndata pickles (476 files,
live legacy path, not staged) instead of reading precomputed data, and
uses a third niche name (`tumor_CAF1(CD105High)` vs. the older source's
`tumor_CAF1_lymphocytes`). A genuine source bug
(`reset_index(inplace=True)` destroying a MultiIndex the rest of the
script needs) was found and fixed (verified logic-preserving). The
long-running compute (2M-row interaction counting, unvectorized, ~1-2s per
niche per sample x 476 samples) died silently partway through (162/476
samples) and was not restarted — status unknown, needs a rerun to
determine if it completes and what it outputs. `figure7def_circos_plots.py`
remains the working, validated source for the published Fig 7d-f panels.

# Issue: myCAF and niche-number identification

`stromal-CAF1(CD105+)` = myCAF, confirmed directly from paper text ("CD105high
are annotated as myCAFs," page 8) — corrects an earlier wrong guess of
`stromal-CAF2(AR+)`.

Niche numbers 1-18 (as used in the paper text/figures) aren't stored in
any data file — inferred from the fixed `niche_order` list used
consistently across `figure5_niche_heatmap.R`/`niche_composition.py`
(position 1 = `luminal`, ..., position 18 = `TLS`). Independently confirmed
for niche 6 (`tumorERG+p53+_ProlifLuminal`, matches Fig 6b) and niches
16-18 (matches Fig 6e's `cols_inflamed` and legend). Niches 2
(`luminal_infiltrated`), 8 (`tumor_CAF1_lymphocytes`), and 9
(`luminal_CAF1(CD105High)`) — used for Fig 7b/S6b-c — are consistent with
their paper-text descriptions but not independently confirmed the same way.

# Issue: don't substitute the look-alike `colormaps.yaml`

`/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/colormaps.yaml`
exists at the top level of the staged data tree and looks like a newer
copy of `resources/colormaps.yaml`, but isn't — confirmed via `diff` it's
missing the `niche`/`metagroups` sections and has different hex values for
shared keys. `resources/colormaps.yaml` (checked into the repo) is the
correct one; don't substitute the other in.

# Issue: `clusters_annotated.parquet` (no `_v2`) is intentionally used by the circos pipeline

`figure7_full_circos_plots.py` reads `clusters_annotated.parquet` (not
`_v2`) — this looks like a stale-data bug at first glance (its `niche`
column disagrees with the verified-correct `_v2` file on 73% of rows), but
isn't: the legacy `compute_interactions.py`/`visualize_interactions_lfc.py`
pipeline genuinely computed the published results using this older
assignment throughout, remapping only output *filenames* to `_v2` niche
names at the very end. Substituting `_v2` into the computation itself
would be the actual deviation from what produced the paper's results.

# Issue: `cell_annotation.parquet` / `niche_frequencies_per_tma_id.parquet` provenance

Both were confirmed genuinely missing (checked `PCa_NHood`, `PCA_NHOODs_clean`,
and their newer staged siblings) until Melissa supplied
`data/melissa-transfer/{clusters_annotated_v2,props_niche_tma_id}.parquet`
and identified they're the same data under a different name — verified
byte-identical via `md5sum` before wiring in (`data/cell_annotation.parquet`
and `LEGACY_DATA_DIR/5-niches/frequencies/niche_frequencies_per_tma_id.parquet`).
This unblocked five KM scripts (6b, 6e, 7b, 7c, S6b-c) and Fig 6a.
