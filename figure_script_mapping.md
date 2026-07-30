# Figure → authoritative old-repo script mapping

**Purpose**: before writing any more `scripts/figures/figureN_*` files, pin
down exactly which script in the old repo
(`/work/FAC/FBM/DBC/mrapsoma/prometex/projects/PCa/`) actually produced each
paper panel. That old repo is the authoritative source — this repo's
`scripts/{02-umaps,03-survival,04-heatmaps,11_niches}/` are a **migration
copy** of it (verified below), not an independent implementation, and the
migration is what introduced the `prepare_data()`/`normalize()` bug and
every hardcoded-path problem documented in `missing_files.md`. Going
forward, new `scripts/figures/figureN_*` files should be written fresh from
the old-repo column, not from the `scripts/{02-umaps,...}` copies.

**Verified**: `000_paper/` in the old repo is byte-identical, file-for-file,
to what's currently sitting in this repo's `scripts/02-umaps/`,
`scripts/03-survival/`, `scripts/04-heatmaps/`, `scripts/11_niches/` (checked
`risk_groups_label.R`, `survival.r`, `2-cell-types-heatmap.R`,
`00_kmeans_clustering.py`, `01_annotation_v2.py` — `diff` clean on all,
`02_umaps/*.py` differ only by two import-path renames). So `000_paper/` is
confirmed as the direct parent of everything already in this repo's
non-`figures/` script directories.

Legend: **Old repo (authoritative)** = path under
`/work/FAC/FBM/DBC/mrapsoma/prometex/projects/PCa/`. **This repo (migrated
copy)** = the already-existing, not-yet-fixed copy under `scripts/` here,
shown only for cross-reference. **New target** = the `scripts/figures/`
file to create (✅ = already exists and already validated against real data
in an earlier branch).

---

## Figure 1 — workflow schematic, panel table, representative images

Not code-derived (BioRender + raw image crops). No script found in the old
repo either. **Gap**, not reconstructable.

## Figure 2 — (a) 34-cell-type heatmap, (b) UMAP of all cells

| Panel | Old repo (authoritative) | This repo (migrated copy) | New target |
|---|---|---|---|
| (a) heatmap | `000_paper/04_heatmaps/2-cell-types-heatmap.R` | `scripts/04-heatmaps/2-cell-types-heatmap.R` | `scripts/figures/figure2_cell_type_heatmap.R` ✅ |
| (b) UMAP | `000_paper/02_umaps/0-umaps.py`, `0-umaps-main-types.py` | `scripts/02-umaps/0-umaps*.py` | `scripts/figures/figure2_umap.py` ✅ |

Note: `000_paper/04_heatmaps/2-1-cell-types-heatmap.R` also exists —
confirmed **older/stale** variant (still on `me3312`'s Mac paths, missing a
later "sample_size" code block present in `2-cell-types-heatmap.R`). Use
`2-cell-types-heatmap.R`, not `2-1-...`.

## Figure 3 — (a) CAF UMAP, (b) CAF subcluster heatmap, (c-e) representative ROIs

| Panel | Old repo (authoritative) | New target |
|---|---|---|
| (a) UMAP | `000_paper/02_umaps/0-umaps-cafs.py` | `scripts/figures/figure3_caf_umap.py` ✅ |
| (b) heatmap | Reuses `2-cell-types-heatmap.R` filtered to CAF labels (already verified while implementing) | `scripts/figures/figure3_caf_heatmap.R` ✅ |
| (c-e) representative ROIs | Not found | **Gap** |

## Figure 4 — (a) patient composition clusters P1–P6, (b) metagroup distribution, (c) KM by patient group, (d-e) Cox PH

| Panel | Old repo (authoritative) | New target |
|---|---|---|
| (a) P1-P6 clustering | **No script found** — the real cluster assignment was read from a precomputed file (`metadata_with_dendrogram_colors_label_pat_id.parquet`, manually labeled), not computed. Confirmed in the `figure-4-composition-survival` branch investigation. | `scripts/figures/figure4_patient_clustering.py` ✅ (reconstructed from Methods text — flagged as an approximation, not a script port) |
| (b) metagroup distribution | Not yet located — check `000_paper/11_niches/111_heatmaps/patient_heatmap.R` / `proportion_heatmap.R` for a metagroup panel, or `100_other_visualization/plot_stacked_frequencies.py` | not yet built |
| (c) KM by patient group | `000_paper/100_other_visualization/risk_groups_label.R` **or** `000_paper/11_niches/113_survival/risk_groups_label.R` — **these are two DIFFERENT scripts, not duplicates** (see note below) | `scripts/figures/figure4_survival.R` ✅ (KM panel currently flagged unresolved due to the P1-P6 cluster-cut mismatch) |
| (d-e) Cox PH | `000_paper/03_survival/survival.r`, `survival-proportions.r`, `plot-cox-hazard-ratio.r`, `survival-cell-freq-groups.R` | `scripts/figures/figure4_survival.R` ✅ (validated, matches paper exactly) |

**Flag**: `000_paper/100_other_visualization/risk_groups_label.R` and
`000_paper/11_niches/113_survival/risk_groups_label.R` have the **same
filename but different content** — the former does generic "risk group
splits" KM analysis (has a `TODO: load needed libraries` header, looks less
finished), the latter does the P1-P6/inflammation/stromogenic/gleason KM
panels and the `clinical.parquet` consistency check. Needs a decision on
which (if either) maps to Fig 4c specifically — flagging rather than
guessing.

## Figure 5 — (a) 18-niche z-scored heatmap, (b) niche correlation matrix, (c) representative cores

| Panel | Old repo (authoritative) | This repo (migrated copy) | New target |
|---|---|---|---|
| clustering (upstream) | `000_paper/11_niches/110_analysis/00_kmeans_clustering.py` (k=24, k-means++) | `scripts/11_niches/110_analysis/00_kmeans_clustering.py` | `scripts/figures/figure5_niche_clustering.py` (not yet created) |
| annotation (upstream) | `000_paper/11_niches/110_analysis/01_annotation_v2.py` | `scripts/11_niches/110_analysis/01_annotation_v2.py` | same file or `figure5_niche_annotation.py` |
| (a) z-score heatmap | `000_paper/11_niches/111_heatmaps/z_score_heatmap.R` | `scripts/11_niches/111_heatmaps/z_score_heatmap.R` | `scripts/figures/figure5_niche_heatmap.R` |
| (b) niche correlation | `000_paper/11_niches/111_heatmaps/niche_pairwise_corrleation.R` | `scripts/11_niches/111_heatmaps/niche_pairwise_corrleation.R` | `scripts/figures/figure5_niche_correlation.R` |
| (c) representative cores | Not found | — | **Gap** |
| Suppl. Fig 7 (50-seed ARI robustness sweep) | Data exists (`.../PCa_NHood/robustness/*_kmeans_robustness.pkl`) but no plotting script found yet | — | **Gap**, may need writing from the raw sweep output |

## Figure 6 — (a) niche abundance heatmap, (b) KM niche 6, (c) composition barplot, (d) niche vs inflammation, (e) immune-niche risk score KM

| Panel | Old repo (authoritative) | New target |
|---|---|---|
| (a) abundance heatmap | `000_paper/11_niches/111_heatmaps/patient_heatmap.R` and/or `proportion_heatmap.R` | `scripts/figures/figure6_niche_abundance_heatmap.R` |
| (b) KM niche 6 | `000_paper/11_niches/113_survival/niche_kaplan_meier_binary.R` | `scripts/figures/figure6_km_niche6.R` |
| (c) composition barplot | `000_paper/11_niches/111_heatmaps/visualize_composition.py` (or `000_paper/100_other_visualization/plot_stacked_frequencies.py` — need to check which one this panel actually is) | `scripts/figures/figure6_composition_barplot.py` |
| (d) niche vs inflammation | `000_paper/11_niches/112_violinplots/inflammation_vis.R` | `scripts/figures/figure6_inflammation_violin.R` |
| (e) immune-niche risk score KM | `000_paper/11_niches/113_survival/risk_groups_label.R` and/or `label_kaplan_meier_binary.R`, `inflammation_outcome.R` | `scripts/figures/figure6_risk_score_km.R` |

## Figure 7 — (a) stromal niche vs stromogenic status, (b-c) KM niche 9 vs myCAF, (d-f) circos interaction plots, (g-h) representative cores

| Panel | Old repo (authoritative) | New target |
|---|---|---|
| (a) stromogenic status | `000_paper/11_niches/112_violinplots/stromogenic_vis.R`, `pairwise_testing_niches.R` | `scripts/figures/figure7_stromogenic_vis.R` |
| (b-c) KM niche 9 vs myCAF | `000_paper/11_niches/113_survival/stromogenic_outcome.R` | `scripts/figures/figure7_km_niche9_mycaf.R` |
| (d-f) circos interaction plots | `000_paper/11_niches/114_interactions/compute_interactions.py`, `circos_plots.py`, `visualize_circos_plot.py`, `nhood_interactions.py`, `visualize_interactions_lfc.py` | `scripts/figures/figure7_circos_interactions.py` |
| (g-h) representative cores | Not found | **Gap** |

---

## Open questions before I start copying/implementing

1. **Fig 4c**: which `risk_groups_label.R` (there are two, different content)?
2. **Fig 6c**: `visualize_composition.py` vs `plot_stacked_frequencies.py` — which one is the actual composition-barplot panel?
3. Should I keep working figure-by-figure (5, then 6, then 7), each its own
   git branch off `figure-4-composition-survival` like before, or something
   else this time given we're starting from a cleaner mapping?
4. For panels marked **Gap** (representative cores, Suppl. Fig 7, Fig 4b
   metagroup distribution not yet located) — reconstruct from Methods text
   like Fig 4a was, or leave flagged and move on?
