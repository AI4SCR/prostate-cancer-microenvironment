"""Hand-maintained table of every in-scope figure/data script for the
regression harness: which paper figure/panel it produces (cross-referenced
against ../../figures.md, the authoritative source), what output files to
hash, and whether it's cheap enough to run by default.

Not parsed out of figures.md automatically -- explicitness here is worth
the duplication; figures.md's `Output` column is the human-readable source
of truth this was built from, kept in sync by hand.

Every entry's `outputs` is a list of glob patterns, relative to the
script's `output_root` ("figures" -> OUTPUT_FIGURES_DIR, "revision" ->
OUTPUT_FIGURES_DIR.parent / "revision"). Globs, not exact filenames,
because several scripts write one file per marker/niche/label and an exact
list would drift every time the underlying data's categories change.

`figure7_full_circos_plots.py` is deliberately absent -- its data was
dropped from DATA_DIR (see data/assets.md), it is not supported.
`scripts/upload-data-to-zenodo.py` is outside this harness's scope (not a
figure/data script).
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Entry:
    script: str  # path relative to repo root
    outputs: list[str]  # glob patterns relative to output_root
    output_root: str = "figures"  # "figures" | "revision"
    slow: bool = False  # excluded from the default fast run
    depends_on: str | None = None  # another entry's script path that must run first


MANIFEST: list[Entry] = [
    # --- Figure 2 ---
    Entry("scripts/figures/figure2_cell_type_heatmap.R", ["figure2/figure2a_cell_type_heatmap.pdf"]),
    Entry("scripts/figures/figure2_umap.py", ["figure2/label=*.pdf", "figure2/value=*.pdf"], slow=True),
    # --- Figure 3 ---
    Entry(
        "scripts/figures/figure3_caf_umap.py",
        ["figure3/excl_markers/*.pdf", "figure3/caf_markers_only/*.pdf", "figure3/stromal/*.pdf"],
        slow=True,
    ),
    Entry("scripts/figures/figure3_caf_heatmap.R", ["figure3/figure3b_caf_heatmap.pdf"]),
    # --- Figure 4 ---
    Entry(
        "scripts/figures/figure4_patient_clustering.py",
        ["figure4/figure4a_stacked_barplot.pdf", "figure4/figure4a_patient_composition.parquet"],
    ),
    Entry(
        "scripts/figures/figure4_metagroup_barplot.py",
        ["figure4/figure4b_metagroup_barplot.pdf"],
        depends_on="scripts/figures/figure4_patient_clustering.py",
    ),
    Entry("scripts/figures/figure4c_patient_cluster_km.R", ["figure4/figure4c_survival_by_patient_cluster.png"]),
    Entry("scripts/figures/figure4de_cox_hazard_ratio.R", ["figure4/figure4_cox_d.png", "figure4/figure4_cox_e.png"]),
    # --- Figure 5 (clustering/annotation chain feeds the heatmap + correlation panels) ---
    Entry(
        "scripts/figures/figure5_niche_clustering.py",
        ["figure5/clusters.parquet", "figure5/*_heatmap.png"],
        slow=True,
    ),
    Entry(
        "scripts/figures/figure5_niche_annotation.py",
        ["figure5/niche_annotations_v2.csv", "figure5/clusters_annotated_v2.parquet"],
        depends_on="scripts/figures/figure5_niche_clustering.py",
        # inherits slow=True: needs clusters.parquet, which only the slow
        # k-means entry above produces -- can't run standalone in the fast subset
        slow=True,
    ),
    Entry("scripts/figures/figure5_niche_heatmap.R", ["figure5/figure5a_niche_zscore_heatmap.pdf"]),
    Entry("scripts/figures/figure5_niche_correlation.R", ["figure5/figure5b_niche_correlation_heatmap.pdf"]),
    # --- Figure 6 ---
    Entry("scripts/figures/figure6_niche_abundance_heatmap.R", ["figure6/figure6a_niche_proportion_heatmap_tma.pdf"]),
    Entry("scripts/figures/figure6_km_niche6.R", ["figure6/niches/*.pdf"]),
    Entry("scripts/figures/figure6c_niche_composition_filtered.py", ["figure6/figure6c_niche_composition_barplot.pdf"]),
    Entry("scripts/figures/figure6_inflammation_violin.R", ["figure6/violin_boxplot_inflammation.pdf"]),
    Entry("scripts/figures/figure6e_immune_risk_score_km.R", ["figure6/kaplan_meier_inflammation_*.pdf"]),
    # --- Figure 7 ---
    Entry(
        "scripts/figures/figure7a_stromogenic_violin.R",
        ["figure7/violin_boxplot_stromogenic.pdf", "figure7/violin_boxplot_stromogenic_sep.pdf"],
    ),
    Entry("scripts/figures/figure7b_niche9_km.R", ["figure7/niches/*luminal_CAF1(CD105High)*.pdf"]),
    Entry("scripts/figures/figure7c_myCAF_km.R", ["figure7/cell_types/*.pdf"]),
    Entry(
        "scripts/figures/figure7def_circos_plots.py",
        [
            "figure7/figure7d_luminal_infiltrated_circos_plot.pdf",
            "figure7/figure7e_tumor_CAF1_lymphocytes_circos_plot.pdf",
            "figure7/figure7f_luminal_CAF1(CD105High)_circos_plot.pdf",
        ],
    ),
    # --- Figure S1 ---
    Entry(
        "scripts/figures/figureS1a_cohort_summary.R",
        ["figureS1/tma-level/*.pdf", "figureS1/patient-level/*.pdf"],
    ),
    Entry("scripts/figures/figureS1b_gleason_concordance.R", ["figureS1/heatmap_gleason_grp_by_gs_grp.pdf"]),
    # --- Figure S2 ---
    Entry(
        "scripts/figures/figureS2_compartment_umap.py",
        ["figureS2/immune/*.pdf", "figureS2/epithelial/*.pdf", "figureS2/endothelial/*.pdf"],
        slow=True,
    ),
    # --- Figure S3 ---
    Entry("scripts/figures/figureS3a_tma_stacked_barplot.py", ["figureS3/figureS3a_stacked_barplot.pdf"]),
    Entry("scripts/figures/figureS3b_cluster_concordance.R", ["figureS3/heatmap_cluster_group_tma_by_patient_cluster_group.pdf"]),
    Entry("scripts/figures/figureS3c_progression_km.R", ["figureS3/progr_patient_cluster_group_final_all.pdf"]),
    # --- Figure S4 ---
    Entry("scripts/figures/figureS4b_niche_mean_composition.py", ["figureS4/figureS4b_niche_composition_barplot.pdf"]),
    Entry("scripts/figures/figureS4b_niche_median_composition.py", ["figureS4/figureS4b_niche_composition_barplot_median.pdf"]),
    # --- Figure S5 ---
    Entry("scripts/figures/figureS5a_p53_niche_barplot.py", ["figureS5/p53_immune_full_barplot_annotated.pdf"]),
    Entry("scripts/figures/figureS5b_inflammation_km.R", ["figureS5/inflammation/*.pdf"]),
    # --- Figure S6 ---
    Entry("scripts/figures/figureS6a_stromogenic_km.R", ["figureS6/stromogenic/*.pdf"]),
    Entry("scripts/figures/figureS6bc_niche_km.R", ["figureS6/niches/*.pdf"]),
    # --- Figure S7 ---
    Entry("scripts/figures/figureS7_ari_robustness.py", ["figureS7/ari_boxplot.png"], slow=True),
    # --- Revision (post-review addenda, not in figures.md's main table) ---
    Entry("scripts/revision/figure6_cox_niche6_erg_p53.R", ["figure6_niche6/*.csv", "figure6_niche6/*.png"], output_root="revision"),
    Entry("scripts/revision/figure6_cox_niche6_gleason.R", ["figure6_niche6/*.csv", "figure6_niche6/*.png"], output_root="revision"),
    Entry("scripts/revision/figure7_cox_niche9_caf1cd105.R", ["figure7_niche9/*.csv", "figure7_niche9/*.png"], output_root="revision"),
    Entry("scripts/revision/figure7_cox_niche9_gleason.R", ["figure7_niche9/*.csv", "figure7_niche9/*.png"], output_root="revision"),
    Entry(
        "scripts/revision/figureS5b_tb_mixed_roi_markers.py",
        ["figureS5b_tb_mixed_roi_markers/*.png", "figureS5b_tb_mixed_roi_markers/*.pdf"],
        output_root="revision",
        slow=True,
    ),
]
