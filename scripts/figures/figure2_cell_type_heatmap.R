# Reproduce Figure 2a: clustered heatmap of normalized mean marker expression
# across the 34 annotated cell types.
#
# Reads the tables scripts/00-data-export/export_for_r.py produces in
# EXPORT_DIR (never BASE_DIR -- see REPRODUCIBILITY.md). Writes to
# EXPORT_DIR/figures/figure2/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

base_dir <- Sys.getenv("BASE_DIR")
export_dir <- Sys.getenv("EXPORT_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(export_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

save_dir <- file.path(export_dir, "figures", "figure2")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# non-biological channels present in intensity_normalized.parquet: DNA
# intercalator + segmentation-kit channels, plus FAP (excluded by the paper's
# Methods due to non-specific staining after in-house conjugation).
non_marker_cols <- c("dna1", "dna2", "icsk1", "icsk2", "icsk3", "fap")
index_cols <- c("sample_id", "object_id", "slide_code", "donor_block_id", "pat_id")

metadata <- read_parquet(file.path(export_dir, "metadata.parquet"))
intensity <- read_parquet(file.path(export_dir, "intensity_normalized.parquet"))

marker_cols <- setdiff(colnames(intensity), c(non_marker_cols, index_cols))
stopifnot("expected 34 markers" = length(marker_cols) == 34)

cells <- intensity |>
  select(sample_id, object_id, all_of(marker_cols)) |>
  inner_join(metadata |> select(sample_id, object_id, label, main_group), by = c("sample_id", "object_id"))

n_cell_types <- n_distinct(cells$label)
stopifnot("expected 34 annotated cell types" = n_cell_types == 34)

mean_expression <- cells |>
  group_by(label) |>
  summarise(across(all_of(marker_cols), mean), main_group = first(main_group), .groups = "drop")

mat <- as.matrix(mean_expression[, marker_cols])
rownames(mat) <- mean_expression$label

main_group_colors <- structure(
  circlize::rand_color(n_distinct(mean_expression$main_group), luminosity = "bright"),
  names = sort(unique(mean_expression$main_group))
)
row_annotation <- rowAnnotation(
  compartment = mean_expression$main_group,
  col = list(compartment = main_group_colors)
)

png(file.path(save_dir, "figure2a_cell_type_heatmap.png"), width = 2400, height = 2600, res = 220)
Heatmap(
  mat,
  name = "mean expr.\n(arcsinh + min-max)",
  right_annotation = row_annotation,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "average",
  clustering_distance_rows = "euclidean",
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = "Figure 2a -- mean marker expression per cell type (n=34)"
)
dev.off()

write_parquet(mean_expression, file.path(save_dir, "figure2a_mean_expression.parquet"))
cat("Saved figure 2a heatmap to", save_dir, "\n")
