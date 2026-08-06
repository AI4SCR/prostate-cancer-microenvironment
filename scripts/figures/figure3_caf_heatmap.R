library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(yaml)

base_dir <- Sys.getenv("BASE_DIR")
export_dir <- Sys.getenv("EXPORT_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(export_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

save_dir <- file.path(output_figures_dir, "figure3")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(export_dir, "intensity_normalized.parquet")
metadata_path <- file.path(export_dir, "metadata.parquet")
colormap_path <- file.path(dirname(export_dir), "resources", "colormaps.yaml")
index_names <- c("sample_id", "object_id")
exclude_channels <- c("fap", "dna1", "dna2", "icsk1", "icsk2", "icsk3")

heatmap.caf <- function(data_path, metadata_path, colormap_path,
                         index_names = c("sample_id", "object_id"),
                         cluster_rows = TRUE) {
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

  filter1 <- grepl("CAF", meta$label)
  # filter2 = grepl('pericytes', meta$label)
  # filter3 = grepl('stromal-', meta$label)
  filter_ <- filter1
  mat <- mat[filter_, ]
  meta <- meta[filter_, ]

  caf_markers <- c(
    "vimentin", "collagen1", "cd146", "cnn1", "smooth_muscle_actin", "cd105", "ar", "pdpn", "egr1", "ces1", "yap1", "beta_catenin"
  )
  filter1 <- colnames(mat) %in% caf_markers
  mat <- mat[, filter1]
  mat <- mat[, caf_markers]

  # AGGREGATION (mean per group)
  aggregate_by <- "label"
  group <- meta[[aggregate_by]]
  mat <- rowsum(mat, group) / as.vector(table(group))
  meta <- meta[!duplicated(group), , drop = FALSE]
  stopifnot(all(meta[[aggregate_by]] == rownames(mat)))

  # Define colour mappings for annotations
  label_colors <- unlist(colormaps$label)
  main_group_colors <- unlist(colormaps$main_group)

  pattern <- paste0("stromal", "-")
  meta$label <- str_replace_all(meta$label, pattern, "")
  names(label_colors) <- str_replace_all(names(label_colors), pattern, "")

  col_anno <- HeatmapAnnotation(
    label = meta$label,
    col = list(label = label_colors),
    na_col = "#F0F0F0",
    show_legend = c(label = TRUE)
  )

  mat <- scale(mat) # z-score each marker (column) across the CAF subclusters
  heatmap_obj <- Heatmap(
    t(mat), # transpose = rotate 90 degrees
    name = "Protein intensity",
    col = circlize::colorRamp2(c(-2, 0, 2), c("lightblue", "white", "lightcoral")),
    cluster_rows = FALSE, # swap
    cluster_columns = TRUE, # swap
    show_row_names = TRUE, # was columns
    show_column_names = FALSE, # was rows
    column_title = NULL, # was row_title
    top_annotation = col_anno, # was right_annotation
    heatmap_legend_param = list(title = "Intensity")
  )

  draw(heatmap_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
}

save_path <- file.path(save_dir, "figure3b_caf_heatmap.pdf")
pdf(save_path, width = 10, height = 10)
heatmap.caf(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  index_names = index_names,
  cluster_rows = TRUE
)
dev.off()
cat("Saved Figure 3b heatmap to", save_path, "\n")
