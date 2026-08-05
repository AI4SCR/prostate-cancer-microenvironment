# Figure → script mapping

Which script produces each panel of each main-text figure. Blank "Script" = not implemented in this repo yet; "—" = not code-derived / no source found. "Validated" = output visually confirmed to match the published panel. "Ported from" = the legacy repo file(s) this script is a 1:1 port of (paths relative to `PCa/`) -- for rows with no `Script` yet, this is our best-guess candidate source (see `open-questions.md` for confidence/caveats), not a confirmed port.

| Figure | Panel | Description | Script | Ported from | Validated |
|---|---|---|---|---|---|
| 1 | a-c | Workflow schematic, antibody panel, representative IMC images | — (BioRender + raw image crops, not code-derived) | — | |
| 2 | a | 34-cell-type heatmap | `scripts/figures/figure2_cell_type_heatmap.R` | `000_paper/04_heatmaps/2-cell-types-heatmap.R` | |
| 2 | b | UMAP of all cells (compartment/cell type/patient/markers) | `scripts/figures/figure2_umap.py` (+ `scripts/data/figure2_umap_embedding.py`) | `000_paper/02_umaps/0-umaps.py` | yes |
| 3 | a | UMAP of CAF cells by subcluster | `scripts/figures/figure3_caf_umap.py` (+ `scripts/data/figure3_caf_umap_embedding.py`) | `000_paper/02_umaps/0-umaps-cafs.py` | |
| 3 | b | CAF subcluster marker-expression heatmap | `scripts/figures/figure3_caf_heatmap.R` | `000_paper/04_heatmaps/2-cell-types-heatmap.R` (`heatmap.caf()`) | |
| 3 | c-e | Representative CAF/fibrocyte ROIs | not implemented | | |
| 4 | a | Patient-level hierarchical clustering (P1-P6) | `scripts/figures/figure4_patient_clustering.py` | `000_paper/sync_paper/05-heterogeneity/stacked-frequencies-label.py` (`group_var=='pat_id'` branch) | yes |
| 4 | b | Mean metagroup distribution per patient cluster | `scripts/figures/figure4_metagroup_barplot.py` | `000_paper/sync_paper/05-heterogeneity/stacked-frequencies-label.py` (tail block) | yes |
| 4 | c | KM survival by patient cluster | `scripts/figures/figure4c_patient_cluster_km.R` | `000_paper/sync_paper/03_survival/patient_risk_group_km.R` | yes |
| 4 | d-e | Cox PH: disease progression / overall survival, per cell type | `scripts/figures/figure4de_cox_hazard_ratio.R` | `000_paper/03_survival/survival-proportions.r` | yes |
| 5 | a | Niche x cell-type z-score heatmap | `scripts/figures/figure5_niche_heatmap.R` (+ `figure5_niche_clustering.py`, `figure5_niche_annotation.py`) | `000_paper/11_niches/111_heatmaps/z_score_heatmap.R`; `000_paper/11_niches/110_analysis/00_kmeans_clustering.py`; `000_paper/11_niches/110_analysis/01_annotation_v2.py` | |
| 5 | b | Niche pairwise Spearman correlation | `scripts/figures/figure5_niche_correlation.R` | `000_paper/11_niches/111_heatmaps/niche_pairwise_corrleation.R` | |
| 5 | c | Representative cores + stacked barplots | not implemented -- no legacy plotting code for per-core composition bars, see `open-questions.md` | — | |
| 6 | a | Per-core niche-abundance heatmap | `scripts/figures/figure6_niche_abundance_heatmap.R` | `000_paper/11_niches/111_heatmaps/patient_heatmap.R` | |
| 6 | b | KM survival by niche 6 abundance | `scripts/figures/figure6_km_niche6.R` | `000_paper/11_niches/113_survival/niche_kaplan_meier_binary.R` | |
| 6 | c | Stacked barplot, niches 6/16/17/18 composition (top panel only) | `scripts/figures/figure6c_niche_composition_filtered.py` | `000_paper/sync_paper/06-spatial-niches/composition/niche_composition.py` (mean-composition-per-niche section, filtered to niches 6/16/17/18) | |
| 6 | d | Niche abundance vs inflammation status | `scripts/figures/figure6_inflammation_violin.R` | `000_paper/11_niches/112_violinplots/inflammation_vis.R` | |
| 6 | e | KM by immune-niche risk score | `scripts/figures/figure6e_immune_risk_score_km.R` | `000_paper/11_niches/113_survival/inflammation_outcome.R` (`cols_inflamed` = niches 16/17/18, `n_inflamed` risk_group 0-3) | |
| 7 | a | Stromal niche abundance vs stromogenic status | `scripts/figures/figure7a_stromogenic_violin.R` | `000_paper/11_niches/112_violinplots/stromogenic_vis.R` | |
| 7 | b | KM by niche 9 abundance | `scripts/figures/figure7b_niche9_km.R` | `000_paper/11_niches/113_survival/niche_kaplan_meier_binary.R` (niche 9 instance) | |
| 7 | c | KM by myCAF abundance | `scripts/figures/figure7c_myCAF_km.R` | `000_paper/11_niches/113_survival/label_kaplan_meier_binary.R` (`stromal-CAF1(CD105+)` instance) | |
| 7 | d-f | Circos plots, cell-cell interactions (niches 2/8/9) | `scripts/figures/figure7def_circos_plots.py` | `000_paper/11_niches/114_interactions/visualize_circos_plot.py` (+ `circos_plots.py`, inlined) | |
| 7 | g-h | Representative cores, niche 8 / niche 9 | not implemented | — (representative IMC images, not code-derived) | |
| S1 | a | EMPaCT cohort description: patient metadata/follow-up summary | `scripts/figures/figureS1a_cohort_summary.R` | `000_paper/01_clinical_metadata/clinical.r` | |
| S1 | b | Core-to-patient Gleason group concordance heatmap | `scripts/figures/figureS1b_gleason_concordance.R` | `000_paper/sync_paper/05-heterogeneity/patient-core-heterogeneity.R` (Gleason-concordance section, `heatmap_gleason_grp_by_gs_grp.pdf`) -- older draft copies of just this half exist at `000_paper/01_clinical_metadata/heterogeneity.r` and `000_paper/100_other_visualization/heterogeneity.R`, see `open-questions.md` | |
| S2 | a-c | Marker expression in compartment-specific UMAPs (epithelial/immune/endothelial), colored by cell type and marker intensity | `scripts/figures/figure2_umap.py` (compartment-subset + intensity loops -- already implemented, no new script) | `000_paper/02_umaps/0-umaps.py` | |
| S3 | a | Core-level (TMA-level) cell-type composition clustering, six groups | `scripts/figures/figureS3a_tma_stacked_barplot.py` | `000_paper/sync_paper/05-heterogeneity/stacked-frequencies-label.py` (`group_var=='tma_id'` branch -- TMA-level sibling of Figure 4a's `pat_id` branch) | |
| S3 | b | Core-to-patient cell-composition cluster concordance heatmap | `scripts/figures/figureS3b_cluster_concordance.R` | `000_paper/sync_paper/05-heterogeneity/patient-core-heterogeneity.R` (cluster-concordance section, `heatmap_cluster_group_tma_by_patient_cluster_group.pdf`) | |
| S3 | c | KM progression-free survival by patient group (P1-P6) | `scripts/figures/figureS3c_progression_km.R` | `000_paper/sync_paper/03_survival/patient_risk_group_km.R` -- same script as Figure 4c; per the paper text ("Fig. 4c, Supplementary Fig 3c" cited together for the same patient groups), this is that script's progression-free plot, which legacy computes but never saves (see `figure4c_patient_cluster_km.R`'s docstring) | |
| S4 | a | Representative niche images | not implemented | — (representative IMC images, not code-derived) | |
| S4 | b | Stacked barplot of cell-type composition per niche (all 18) | `scripts/figures/figureS4b_niche_mean_composition.py` | `000_paper/sync_paper/06-spatial-niches/composition/niche_composition.py` (mean-composition-per-niche section) | |
| S5 | a | Stacked barplot of niche 6 cell-type composition, ordered by per-core niche 6 proportion | not implemented -- no legacy plotting code for per-core composition bars, see `open-questions.md` | — | |
| S5 | b | KM overall survival by inflammation status (patient-level clinical variable) | `scripts/figures/figureS5b_inflammation_km.R` | `000_paper/11_niches/113_survival/risk_groups_label.R` (inflammation section, `p_survival_inflam`) | |
| S6 | a | KM progression-free survival by stromogenic status | `scripts/figures/figureS6a_stromogenic_km.R` | `000_paper/11_niches/113_survival/risk_groups_label.R` (stromogenic section, `p_prog_stromo`, `progr_stromogenic.pdf`) | |
| S6 | b-c | KM progression-free survival by niche 2 / niche 8 abundance | `scripts/figures/figureS6bc_niche_km.R` | `000_paper/11_niches/113_survival/niche_kaplan_meier_binary.R` (niche 2 / niche 8 instances) | |
| S7 | a | 50-seed k-means ARI robustness sweep (boxplots of pairwise adjusted Rand index per run; best-agreement run highlighted) | not implemented | — no plotting script found anywhere in the legacy repo (data exists at `.../PCa_NHood/robustness/*_kmeans_robustness.pkl`); would need to be written from raw sweep output | |
