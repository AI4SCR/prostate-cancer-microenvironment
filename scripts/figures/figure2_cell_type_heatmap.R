library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(yaml)
# viridis isn't installed in this environment; viridisLite::inferno() is the
# same function viridis re-exports (viridis depends on viridisLite for it),
# so this is not a substitution of behavior, just of which package name is
# loaded to reach the identical function.
library(viridisLite)

base_dir <- Sys.getenv("BASE_DIR")
export_dir <- Sys.getenv("EXPORT_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(export_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

save_dir <- file.path(output_figures_dir, "figure2")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(export_dir, "intensity_normalized.parquet")
metadata_path <- file.path(export_dir, "metadata.parquet")
colormap_path <- file.path(dirname(export_dir), "resources", "colormaps.yaml")
index_names <- c("sample_id", "object_id")
exclude_channels <- c("fap", "dna1", "dna2", "icsk1", "icsk2", "icsk3")

heatmap.agg <- function(data_path, metadata_path, colormap_path,
                         index_names = c("sample_id", "object_id"),
                         cluster_rows = FALSE, aggregate_by = NULL, column_split = TRUE) {
  # LOAD DATA
  data <- read_parquet(data_path)
  data <- data |> select(-all_of(exclude_channels))
  meta <- read_parquet(metadata_path)
  colormaps <- yaml::read_yaml(colormap_path)

  # SORT
  ord <- order(meta$label)
  data <- data[ord, ]
  meta <- meta[ord, ]

  # INDEX
  data.index <- data[, index_names]
  data <- data[, !(names(data) %in% index_names)]
  meta.index <- meta[, index_names]
  meta <- meta[, !(names(meta) %in% index_names)]
  stopifnot(all(meta.index == data.index))
  index <- do.call(paste0, meta.index[index_names])

  mat <- as.matrix(data)
  rownames(mat) <- index

  # FILTER
  remove <- meta$label == "mix-vessels-PMN-MDSCs"
  mat <- mat[!remove, ]
  meta <- meta[!remove, ]

  # AGGREGATION (mean per group)
  if (!is.null(aggregate_by)) {
    group <- meta[[aggregate_by]]
    mat <- rowsum(mat, group) / as.vector(table(group))
    meta <- meta[!duplicated(group), , drop = FALSE]
    stopifnot(all(meta[[aggregate_by]] == rownames(mat)))
  }

  label_colors <- unlist(colormaps$label)
  main_group_colors <- unlist(colormaps$main_group)

  col_anno <- HeatmapAnnotation(
    label = meta$label,
    main_group = meta$main_group,
    col = list(label = label_colors, main_group = main_group_colors),
    na_col = "#F0F0F0",
    show_legend = c(label = TRUE, main_group = TRUE)
  )

  if (column_split) {
    heatmap_obj <- Heatmap(
      t(mat),
      name = "Protein intensity",
      col = inferno(256),
      cluster_rows = TRUE,
      cluster_columns = cluster_rows,
      show_row_names = TRUE,
      show_column_names = FALSE,
      column_split = meta$main_group,
      column_title = NULL,
      top_annotation = col_anno,
      heatmap_legend_param = list(title = "Intensity")
    )
  } else {
    heatmap_obj <- Heatmap(
      t(mat),
      name = "Protein intensity",
      col = inferno(256),
      cluster_rows = TRUE,
      cluster_columns = cluster_rows,
      show_row_names = TRUE,
      show_column_names = FALSE,
      column_title = NULL,
      top_annotation = col_anno,
      heatmap_legend_param = list(title = "Intensity")
    )
  }

  draw(heatmap_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
}

save_path <- file.path(save_dir, "figure2a_cell_type_heatmap.pdf")
pdf(save_path, width = 15, height = 10)
heatmap.agg(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  index_names = index_names,
  cluster_rows = TRUE,
  aggregate_by = "label",
  column_split = FALSE
)
dev.off()
cat("Saved Figure 2a heatmap to", save_path, "\n")
