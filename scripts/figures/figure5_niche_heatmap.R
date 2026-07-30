# Reproduce Figure 5a: z-scored niche x cell-type composition heatmap.
#
# 1:1 port of the old repo's 000_paper/11_niches/111_heatmaps/z_score_heatmap.R
# (paths only changed; see figure_script_mapping.md). All inputs here are
# precomputed composition/annotation tables with no reproducing script in
# this repo (figure5_niche_annotation.py is blocked on the missing
# niche_annotations_revised.xlsx -- see that script's docstring), so they're
# read from LEGACY_DATA_DIR rather than regenerated.
#
# Writes to $EXPORT_DIR/figures/figure5/.

library(dotenv)
load_dot_env()

library(arrow)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)

export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

save_dir <- file.path(export_dir, "figures", "figure5")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

niches_dir <- file.path(legacy_dir, "5-niches")
composition_dir <- file.path(niches_dir, "visualization", "composition")

df_annotated <- read_parquet(file.path(niches_dir, "annotation", "clusters_annotated_v2.parquet"))
df_zscore <- read_parquet(file.path(composition_dir, "niche_heatmap_data.parquet"))
df_means <- read_parquet(file.path(composition_dir, "mean_celltype_composition_per_niche.parquet"))
df_medians <- read_parquet(file.path(composition_dir, "median_celltype_composition_per_niche.parquet"))

matrix <- df_zscore %>%
  column_to_rownames("niche") %>%
  as.matrix()

info_niches <- read.csv(file.path(niches_dir, "annotation", "niche_annotations_v2.csv"))

info_niches <- info_niches %>%
  select(-cluster) %>%
  distinct() %>%
  filter(niche %in% rownames(matrix)) %>%
  filter(niche != "unassigned")

niche_order <- c(
  "luminal", "luminal_infiltrated", "basal_luminal_glands",
  "tumor(ERG+)", "tumor(ERG+)_luminal", "tumorERG+p53+_ProlifLuminal", "tumor_CAF1(CD105High)", "tumor_CAF1_lymphocytes",
  "luminal_CAF1(CD105High)", "luminal_CAF1(CD105-)", "bloodvessels_CAF1(CD105-)",
  "CAF1_CD105-_infiltrated", "CAFs_lymphocytes", "CAF2s_enriched", "CAF2(AR-)_enriched", "immune_bloodvessels_CAF1(CD105-)",
  "Macrophages_Tcells_CAF1(CD105-)", "TLS"
)
order <- match(niche_order, info_niches$niche)
info_niches <- info_niches[order, ]

info_celltypes <- df_annotated %>%
  select(label, main_group) %>%
  distinct() %>%
  filter(main_group != "undefined") %>%
  arrange(label)

matrix <- matrix[info_niches$niche, ]
matrix <- matrix[, info_celltypes$label]

means <- df_means %>%
  column_to_rownames("niche") %>%
  as.matrix()
means <- means[info_niches$niche, info_celltypes$label]

# --- annotation colors ---
niche_colors <- setNames(unique(info_niches$niche_color), unique(info_niches$niche))
meta_niche_colors <- setNames(unique(info_niches$meta_niche_color), unique(info_niches$meta_niche))

row_ha <- rowAnnotation(
  df = data.frame(
    meta_niche = factor(info_niches$meta_niche, levels = names(meta_niche_colors)),
    niche = factor(info_niches$niche, levels = names(niche_colors))
  ),
  col = list(
    meta_niche = meta_niche_colors,
    niche = niche_colors
  ),
  show_legend = c(FALSE, FALSE)
)

cell_colors <- setNames(
  c("#8c510a", "#d73027", "#1a9850", "#4575b4", "#bdbdbd"),
  c("endothelial", "epithelial", "immune", "stromal", "undefined")
)

celltype_colors <- setNames(
  c(
    "#05acc6", "#a877ac", "#ff28fd", "#bc9157", "#56642a", "#97ff00", "#ff3464",
    "#a0e491", "#d3008c", "#77c6ba", "#90318e", "#8c9ab1", "#f4bfb1", "#9ee2ff",
    "#829026", "#b8ba01", "#d60000", "#ff8ec8", "#366962", "#ff6200", "#009e7c",
    "#9a6900", "#00c846", "#ae083f", "#c86e66", "#00fdcf", "#79525e", "#ffa52f",
    "#afa5ff", "#018700", "#f2cdff", "#b500ff", "#953f1f", "#fdf490", "#bdbdbd"
  ),
  c(
    "endothelial-blood-vessel(ERG+)", "endothelial-blood-vessel(ERG-)", "endothelial-lymphatic",
    "epithelial-(ERG+CD44+)", "epithelial-basal", "epithelial-luminal", "epithelial-luminal(ERG+)",
    "epithelial-luminal(ERG+p53+)", "epithelial-luminal(Ki67+)", "epithelial-luminal(p53+)", "epithelial-neuroendocrine",
    "epithelial-transient", "immune-B-cells(CD20+)", "immune-BM-derived-fibrocytes",
    "immune-PMN-MDSCs(CD11b+CD66b+)", "immune-T-cells(CD3+)", "immune-T-cells_cytotoxic(CD3+CD8a+)",
    "immune-T-cells_helper(CD3+CD4+)", "immune-T-cells_regulatory(CD3+CD4+FoxP3+)", "immune-T-helper-B-cells",
    "immune-T-helper-T-cytotoxic", "immune-T-helper-macrophages", "immune-macrophages(CD68+)", "mix-vessels-PMN-MDSCs",
    "stromal-(Ki67+)", "stromal-CAF1(CD105+)", "stromal-CAF1(CD105-)", "stromal-CAF1(CD105-EGR1+)",
    "stromal-CAF2(AR+)", "stromal-CAF2(AR+CES1+)", "stromal-CAF2(AR+EGR1+)", "stromal-CAF2(AR-)",
    "stromal-mesenchymal-neuroendocrine", "stromal-pericytes", "undefined"
  )
)

col_ha <- HeatmapAnnotation(
  celltype = info_celltypes$label,
  main_group = info_celltypes$main_group,
  col = list(
    main_group = cell_colors,
    celltype = celltype_colors
  ),
  annotation_name_side = "left"
)

# --- stats annotation ---
df_stats <- read_parquet(file.path(composition_dir, "niche_abundance_stats.parquet"))
df_stats <- df_stats %>% filter(niche != "unassigned")

medians <- df_stats$median_frequency
names(medians) <- df_stats$niche
medians <- medians[rownames(matrix)]

presence <- df_stats$num_samples
names(presence) <- df_stats$niche
presence <- presence[rownames(matrix)]
col_fun_presence <- colorRamp2(range(presence, na.rm = TRUE), c("lightyellow", "darkred"))
my_colors <- col_fun_presence(presence)

ha <- rowAnnotation(
  frequency = anno_numeric(
    medians,
    rg = c(0, max(medians)),
    x_convert = function(x) round(x, 2),
    labels_format = function(x) sprintf("%.2f", x),
    bg_gp = gpar(fill = my_colors, col = "black")
  ),
  annotation_name_rot = 0
)

meta_niche_levels <- c("luminal_other_epithelial_cells", "tumor_cells_+/-_CAFs", "luminal_CAFs", "CAF_subtypes_+/-_immune", "immune_cells")
info_niches$meta_niche <- factor(info_niches$meta_niche, levels = meta_niche_levels)

# --- draw heatmap ---
col_fun <- colorRamp2(c(-3, 0, 3), c("#2166ac", "white", "#b2182b"))
h <- Heatmap(
  matrix,
  col = col_fun,
  name = "Z-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "spearman",
  left_annotation = row_ha,
  right_annotation = ha,
  top_annotation = col_ha,
  column_split = info_celltypes$main_group,
  heatmap_legend_param = list(title = "z-score"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_row_names = TRUE,
  show_column_names = TRUE
)

lgd_presence <- Legend(
  title = "N samples",
  col_fun = col_fun_presence,
  at = pretty(range(presence, na.rm = TRUE), n = 5),
  labels = pretty(range(presence, na.rm = TRUE), n = 5)
)

plot_path <- file.path(save_dir, "figure5a_niche_zscore_heatmap.pdf")
pdf(plot_path, width = 18, height = 14)
draw(h, heatmap_legend_side = "right", annotation_legend_list = list(lgd_presence), annotation_legend_side = "right")
dev.off()
cat("Saved Figure 5a to", plot_path, "\n")
