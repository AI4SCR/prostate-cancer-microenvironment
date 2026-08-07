library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(yaml)

base_dir <- Sys.getenv("BASE_DIR")
data_dir <- Sys.getenv("DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(data_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(data_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

save_dir <- file.path(output_figures_dir, "figure3")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(data_dir, "cells", "intensity_normalized.parquet")
metadata_path <- file.path(data_dir, "cells", "metadata.parquet")
colormap_path <- file.path(dirname(data_dir), "resources", "colormaps.yaml")
index_names <- c("sample_id", "object_id")
exclude_channels <- c("fap", "dna1", "dna2", "icsk1", "icsk2", "icsk3")

heatmap.caf <- function(data_path,
                         metadata_path,
                         colormap_path,
                         index_names = c("sample_id", "object_id"),
                         cluster_rows = TRUE) {
  # -----------------------------------------------------------------------------
  # LOAD DATA
  data <- read_parquet(data_path)
  data <- data |> select(-all_of(exclude_channels))

  meta <- read_parquet(metadata_path)
  colormaps <- yaml::read_yaml(colormap_path)

  # -----------------------------------------------------------------------------
  # SORT
  ord <- order(meta$label)


  data <- data[ord, ]
  meta <- meta[ord, ]

  # -----------------------------------------------------------------------------
  # INDEX
  data.index <- data[, index_names]
  data <- data[, !(names(data) %in% index_names)]

  # meta_stroma = meta %>% filter(grepl("stromal", label))
  ## for each sample calculate total number of stromal cells
  stroma_size <- meta %>%
    filter(grepl("stromal", label)) %>%
    group_by(sample_id) %>%
    summarise(n_stroma = n())

  sample_size <- meta %>%
    group_by(sample_id) %>%
    summarise(n_total = n())

  filter1 <- grepl("CAF", meta$label)

  ## number of each caf type CAFs per sample
  caf_counts <- meta %>%
    filter(filter1) %>%
    group_by(sample_id, label) %>%
    summarise(n_caf = n()) %>%
    ungroup() %>%
    left_join(stroma_size, by = "sample_id") %>%
    mutate(prop_caf = n_caf / n_stroma)

  # caf_counts = meta %>%
  #   filter(filter1) %>%
  #   group_by(sample_id, label) %>%
  #   summarise(n_caf = n()) %>%
  #   ungroup() %>%
  #   left_join(sample_size, by = 'sample_id') %>%
  #   mutate(prop_caf = n_caf / n_total)

  ## make to wide
  caf_wide <- caf_counts %>%
    select(sample_id, label, prop_caf) %>%
    tidyr::pivot_wider(names_from = label, values_from = prop_caf, values_fill = 0)

  median_freq <- apply(caf_wide[, -1], 2, median)
  std_freq <- apply(caf_wide[, -1], 2, sd)
  ## count how many samples have bigger value than median freq
  # freq_above_median = apply(caf_wide[,-1], 2, function(x) sum(x > median_freq[names(x)]))


  meta.index <- meta[, index_names]
  meta <- meta[, !(names(meta) %in% index_names)]

  stopifnot(all(meta.index == data.index))
  index <- do.call(paste0, meta.index[index_names])

  # order = match(df1, rownames(df2))

  mat <- as.matrix(data)
  rownames(mat) <- index

  # -----------------------------------------------------------------------------
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

  # DISCLOSED DEVIATION from `2-1-cell-types-heatmap.R`'s verbatim list (see
  # figures.md's Fig 3b Issue note): matched directly to the marker order
  # visible in the published panel instead -- vimentin, collagen1, cd146,
  # cnn1, smooth_muscle_actin (alpha-SMA), cd105, ar, pdpn, egr1, ces1, yap1
  # (pYAP), beta_catenin. `pdpn` included, `c_casp3`/`ki_67` dropped -- the
  # opposite of the verbatim source's list.
  caf_markers <- c(
    "vimentin",
    "collagen1",
    "cd146",
    "cnn1",
    "smooth_muscle_actin",
    "cd105",
    "ar",
    "pdpn",
    "egr1",
    "ces1",
    "yap1",
    "beta_catenin"
  )
  filter1 <- colnames(mat) %in% caf_markers
  mat <- mat[, filter1]
  mat <- mat[, caf_markers]

  # -----------------------------------------------------------------------------
  # AGGREGATION (mean per group)
  aggregate_by <- "label"
  group <- meta[[aggregate_by]]
  mat <- rowsum(mat, group) / as.vector(table(group))
  meta <- meta[!duplicated(group), , drop = FALSE]
  stopifnot(all(meta[[aggregate_by]] == rownames(mat)))

  # Define colour mappings for annotations
  label_colors <- unlist(colormaps$label)
  main_group_colors <- unlist(colormaps$main_group)
  patient_colors <- unlist(colormaps$pat_id)

  pattern <- paste0("stromal", "-")
  meta$label <- str_replace_all(meta$label, pattern, "")
  names(label_colors) <- str_replace_all(names(label_colors), pattern, "")

  # -----------------------------------------------------------------------------
  # Create row annotations
  cell_order <- c(
    "stromal-CAF1(CD105-)", "stromal-CAF1(CD105-EGR1+)", "stromal-CAF1(CD105+)", "stromal-CAF2(AR+)",
    "stromal-CAF2(AR-)", "stromal-CAF2(AR+CES1+)", "stromal-CAF2(AR+EGR1+)"
  )
  cell_order_meta <- sub("^stromal-", "", cell_order)

  meta <- meta[meta$label %in% cell_order_meta, , drop = FALSE]
  meta$label <- factor(meta$label, levels = cell_order_meta)
  meta <- meta[order(meta$label), , drop = FALSE]
  meta$label <- as.character(meta$label)
  mat <- mat[cell_order, ]

  col_anno <- HeatmapAnnotation(
    label = meta$label,
    col = list(label = label_colors),
    na_col = "#F0F0F0",
    show_legend = c(label = TRUE)
  )

  # # reorder columns to match heatmap order = meta$label
  # caf_box_mat <- caf_wide[, paste("stromal-", meta$label, sep = ""), drop = FALSE]
  # ##remove stromal from colnames
  # colnames(caf_box_mat) <- meta$label
  #
  # bottom_anno <- HeatmapAnnotation(
  #   caf_box = anno_boxplot(
  #     as.matrix(caf_box_mat),
  #     height = unit(4, "cm"),
  #     gp = gpar(fill = "grey60")
  #   ),
  #   annotation_name_side = "left",
  #   which = "column"
  # )
  # cell_order <- c("stromal-CAF1(CD105-)", "stromal-CAF1(CD105-EGR1+)","stromal-CAF1(CD105+)","stromal-CAF2(AR+)",
  #                 "stromal-CAF2(AR-)" ,"stromal-CAF2(AR+CES1+)","stromal-CAF2(AR+EGR1+)")
  # mat <- mat[cell_order, ]
  # -----------------------------------------------------------------------------
  # Draw the heatmap
  # mat = scale(mat)
  mat <- scale(mat)
  heatmap_obj <- Heatmap(
    t(mat), # transpose = rotate 90 degrees
    name = "Protein intensity",
    # col = inferno(256),
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
    cluster_rows = FALSE, # swap
    cluster_columns = FALSE, # swap
    show_row_names = TRUE, # was columns
    show_column_names = FALSE, # was rows
    # column_split = meta$main_group,   # was row_split
    column_title = NULL, # was row_title
    top_annotation = col_anno, # was right_annotation
    # bottom_annotation = bottom_anno,
    heatmap_legend_param = list(title = "Intensity")
  )

  draw(heatmap_obj,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
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
