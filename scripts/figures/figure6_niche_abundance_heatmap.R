library(dotenv)
load_dot_env()

library(arrow)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(entropy)
library(yaml)

export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

save_dir <- file.path(output_figures_dir, "figure6")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
resources_dir <- file.path(dirname(dirname(output_figures_dir)), "resources")
colormaps_path <- file.path(resources_dir, "colormaps.yaml")

freq_path <- file.path(legacy_dir, "5-niches", "frequencies", "niche_frequencies_per_tma_id.parquet")
df_props <- read_parquet(freq_path)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()

matrix <- df_props %>%
  column_to_rownames("tma_id") %>%
  as.matrix()

clinical <- clinical %>%
  filter(tma_id %in% rownames(matrix)) %>%
  arrange(tma_id)

# clusters <- clusters %>%
#   select(tma_id, cluster) %>%
#   distinct()

df_metadata <- clinical %>%
  select(
    pat_id,
    os_status,
    disease_progr,
    # psa_progr,
    # recurrence,
    # gs_grp,
    gleason_grp,
    inflammation,
    stromogenic_smc_loss_reactive_stroma_present,
    # d_amico_risk,
    # glandular_atrophy_pin,
    # ln_status,
    # DISCLOSED DEVIATION from the verbatim source: tma_id is commented out
    # here in the legacy script, but matrix indexing below (`df_metadata$tma_id`)
    # requires it -- a genuine bug in the source itself. Re-enabled.
    tma_id
  ) %>%
  distinct()

# df_metadata <- merge(df_metadata, clusters, by = "tma_id", all.x = TRUE)

matrix <- matrix[as.character(df_metadata$tma_id), ] # reorder matrix rows to match metadata

# hue_list
hue_order_list <- list(
  cause_of_death = c("alive", "PCa_death", "non-PCa_death"),
  os_status = c("alive", "dead"),
  disease_progr = c(0, 1),
  clinical_progr = c(0, 1),
  psa_progr = c(0, 1),
  recurrence = c("no_recurrence", "recurrence"),
  gs_grp = c("1", "2", "3", "4", "5", "nan"),
  gleason_grp = c("1", "2", "3", "4", "5"),
  inflammation = c("yes", "no", "nan"),
  stromogenic_smc_loss_reactive_stroma_present = c("yes", "no", "nan"),
  glandular_atrophy_pin = c("yes", "no", "nan"),
  ln_status = c(0, 1)
)

# Define annotation colors (auto-detect categorical columns)
for (col in names(hue_order_list)) {
  if (col %in% colnames(df_metadata)) {
    df_metadata[[col]] <- factor(df_metadata[[col]], levels = hue_order_list[[col]])
  }
}

## import yaml for colors
custom_annotation_colors <- yaml::read_yaml(colormaps_path)

df_metadata$pat_id <- as.character(df_metadata$pat_id)
annotation_colors <- list()
for (col in colnames(df_metadata)) {
  if (col %in% names(custom_annotation_colors)) {
    # DISCLOSED DEVIATION: yaml::read_yaml() returns a nested list, not the
    # named vector ComplexHeatmap requires; it also parses bareword yes/no
    # keys as YAML 1.1 booleans, not the strings "yes"/"no" the data uses
    # (same gotcha fixed elsewhere in this repo), and colormaps.yaml's
    # gleason_grp keys ("1.0".."5.0") don't match the factor levels
    # produced by factor(<numeric>, ...) ("1".."5"). unlist() + remap.
    vals <- unlist(custom_annotation_colors[[col]])
    names(vals)[names(vals) == "TRUE"] <- "yes"
    names(vals)[names(vals) == "FALSE"] <- "no"
    names(vals) <- sub("\\.0$", "", names(vals))
    annotation_colors[[col]] <- vals
  } else if (col %in% names(hue_order_list)) {
    # Generate consistent colors based on predefined order
    unique_vals <- hue_order_list[[col]]
    annotation_colors[[col]] <- structure(
      circlize::rand_color(length(unique_vals)),
      names = unique_vals
    )
  }
}

# Create HeatmapAnnotation
# tma_id is kept in df_metadata above only for row-indexing (matrix reorder,
# distinct()); the published panel's legend has no tma_id track, so it's
# excluded from the displayed annotation here.
df_metadata_display <- df_metadata %>% select(-tma_id)
annotation_colors_display <- annotation_colors[names(annotation_colors) != "tma_id"]
ha <- rowAnnotation(df = df_metadata_display, col = annotation_colors_display, show_legend = c(FALSE, rep(TRUE, length(df_metadata_display) - 1)))

col_fun <- colorRamp2(c(0, 1), c("white", "darkgreen"))

jsd_distance <- function(x, y) {
  entropy_x <- entropy(x)
  entropy_y <- entropy(y)
  m <- (x + y) / 2

  ### JSD
  jsd <- entropy(m) - (entropy_x + entropy_y) / 2
  return(jsd)
}

# --- define annotation colors ---
# read colormap
colormap_niche <- yaml::read_yaml(colormaps_path)$niche

# ensure named character vector
niche_colors <- unlist(colormap_niche)

# DISCLOSED DEVIATION from the verbatim source: reorder/filter columns via
# niche_annotations_v2.csv (dropping "unassigned"), matching the pattern
# already used in figure5_niche_correlation.R and an earlier version of
# this same legacy script (commit a7d0b60) -- the current sync_paper
# version regressed to sourcing niche colors from colormaps.yaml directly,
# which still lists "unassigned" and reintroduces it as a spurious 19th
# column not present in the published panel.
info_niches <- read.csv(file.path(legacy_dir, "5-niches", "annotation", "niche_annotations_v2.csv"))
info_niches <- info_niches %>%
  select(-cluster) %>%
  distinct() %>%
  filter(niche != "unassigned")

matrix <- matrix[, info_niches$niche]

# annotation dataframe (aligned with matrix columns)
niche_anno_df <- data.frame(
  niche = factor(colnames(matrix), levels = names(niche_colors))
)

# build annotation
niche_anno <- HeatmapAnnotation(
  df = niche_anno_df,
  col = list(niche = niche_colors),
  show_legend = TRUE
)

k <- 5
# Plot heatmap
p <- ComplexHeatmap::Heatmap(matrix,
  col = col_fun,
  name = "Proportion",
  cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_distance_columns = "pearson",
  clustering_distance_rows = jsd_distance,
  show_row_names = FALSE, show_column_names = TRUE,
  right_annotation = ha,
  top_annotation = niche_anno,
  row_split = k,
  row_dend_width = unit(3, "cm"),
  heatmap_legend_param = list(title = "Proportion"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8)
)

plot_path <- file.path(save_dir, "figure6a_niche_proportion_heatmap_tma.pdf")
pdf(plot_path, width = 18, height = 15)
ht <- draw(p)
dev.off()
cat("Saved Figure 6a to", plot_path, "\n")
