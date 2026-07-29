# Reproduce Figure 3b: heatmap of mean marker expression per CAF subcluster.
#
# Paper: "Values represent average z-scored marker expression per cluster."
# Unlike Figure 2a (raw normalized mean, no z-scoring), here each marker
# column is z-scored across the 10 stromal/CAF subclusters after averaging.
#
# Reads the tables scripts/00-data-export/export_for_r.py produces in
# EXPORT_DIR (never BASE_DIR -- see REPRODUCIBILITY.md). Writes to
# EXPORT_DIR/figures/figure3/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(ComplexHeatmap)

base_dir <- Sys.getenv("BASE_DIR")
export_dir <- Sys.getenv("EXPORT_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(export_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

save_dir <- file.path(export_dir, "figures", "figure3")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

non_marker_cols <- readLines(file.path(export_dir, "non_marker_channels.txt"))
index_cols <- c("sample_id", "object_id", "slide_code", "donor_block_id", "pat_id")

metadata <- read_parquet(file.path(export_dir, "metadata.parquet"))
intensity <- read_parquet(file.path(export_dir, "intensity_normalized.parquet"))

marker_cols <- setdiff(colnames(intensity), c(non_marker_cols, index_cols))
stopifnot("expected 34 markers" = length(marker_cols) == 34)

cells <- intensity |>
  select(sample_id, object_id, all_of(marker_cols)) |>
  inner_join(metadata |> select(sample_id, object_id, label, main_group), by = c("sample_id", "object_id")) |>
  filter(main_group == "stromal")

n_caf_labels <- n_distinct(cells$label)
stopifnot("expected 10 stromal/CAF labels" = n_caf_labels == 10)

mean_expression <- cells |>
  group_by(label) |>
  summarise(across(all_of(marker_cols), mean), .groups = "drop")

mat <- as.matrix(mean_expression[, marker_cols])
rownames(mat) <- mean_expression$label
mat_z <- scale(mat)  # z-score each marker (column) across the 10 subclusters

png(file.path(save_dir, "figure3b_caf_heatmap.png"), width = 2200, height = 1400, res = 220, type = "cairo")
Heatmap(
  mat_z,
  name = "z-scored\nmean expr.",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "average",
  clustering_distance_rows = "euclidean",
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "Figure 3b -- z-scored mean marker expression per CAF subcluster (n=10)"
)
dev.off()

write_parquet(mean_expression, file.path(save_dir, "figure3b_mean_expression.parquet"))
cat("Saved figure 3b heatmap to", save_dir, "\n")
