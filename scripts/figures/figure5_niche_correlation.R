library(dotenv)
load_dot_env()

library(arrow)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)

data_dir <- Sys.getenv("DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(data_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))

save_dir <- file.path(output_figures_dir, "figure5")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

niches_dir <- file.path(data_dir, "niches")

df_props <- read_parquet(file.path(niches_dir, "niche_frequencies_per_tma_id_stacked.parquet"))
clinical <- read_parquet(file.path(data_dir, "clinical.parquet"))

matrix <- df_props %>%
  column_to_rownames("tma_id") %>%
  as.matrix()

clinical <- clinical %>%
  filter(tma_id %in% rownames(matrix)) %>%
  arrange(tma_id)

df_metadata <- clinical %>%
  select(
    pat_id, os_status, disease_progr, gs_grp, gleason_grp,
    inflammation, stromogenic_smc_loss_reactive_stroma_present, tma_id
  ) %>%
  distinct()

matrix <- matrix[as.character(df_metadata$tma_id), ]

info_niches <- read.csv(file.path(niches_dir, "niche_annotations.csv"))
info_niches <- info_niches %>%
  select(-cluster) %>%
  distinct() %>%
  filter(niche != "unassigned")

matrix <- matrix[, info_niches$niche]

corr_matrix <- cor(matrix, use = "pairwise.complete.obs", method = "spearman")

col_fun <- colorRamp2(c(-1, 0, 1), c("#2166ac", "white", "#b2182b"))
h <- Heatmap(
  corr_matrix,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(sprintf("%.2f", corr_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 8))
  }
)

plot_path <- file.path(save_dir, "figure5b_niche_correlation_heatmap.pdf")
pdf(plot_path, width = 18, height = 14)
draw(h, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Saved Figure 5b to", plot_path, "\n")
