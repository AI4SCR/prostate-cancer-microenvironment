# Reproduce Figure 6a: per-core (TMA) niche-proportion heatmap.
#
# 1:1 port of the old repo's
# 000_paper/sync_paper/06-spatial-niches/abundance/heatmap_frequencies.R --
# REPLACES an earlier, wrong-script version of this file that ported
# 000_paper/11_niches/111_heatmaps/patient_heatmap.R (patient-level: rows =
# patients, row annotations cause_of_death/clinical_progr/psa_progr/
# recurrence/gs_grp). That was the wrong legacy source: the paper's actual
# Fig 6a legend is "Clustered heatmap of niche abundance scores per tumor
# core" (TMA/core-level), and its row annotations
# (pat_id/os_status/disease_progr/gleason_grp/inflammation/
# stromogenic_smc_loss_reactive_stroma_present) match this script's
# `df_metadata` select exactly, confirmed by direct user identification --
# see open-questions.md. `patient_heatmap.R`'s own separate dendrogram-split
# risk-group KM tail block does not exist in this script and is dropped
# (it never corresponded to any published panel anyway, just a byproduct of
# the wrong-source script).
#
# BLOCKED: this script's required input, `niche_frequencies_per_tma_id.parquet`,
# does not exist anywhere accessible in this environment -- checked
# `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood` (nothing matching
# `niche_frequencies*`) and `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/frequencies/`
# (has a similarly-named `niche_frequencies_per_tma.parquet`, no `_id`, but
# it uses an older pre-"revised annotation" niche naming scheme, e.g.
# `TLS_Bcells_Tcells`, `canonical_BLepithelium` -- not the same data, not a
# usable substitute). Per CLAUDE.md's hard constraint, we do not write a
# script to compute this data ourselves -- this script is ported and
# correct, but will fail at the `read_parquet()` call below until that file
# turns up staged somewhere. See data/assets.md and open-questions.md.
#
# Legacy bug, preserved verbatim per explicit user instruction (no fix,
# strict 1:1 port): `tma_id` is commented out of `df_metadata`'s `select()`
# but still referenced two lines later (`matrix[as.character(df_metadata$tma_id), ]`).
# Since `df_metadata$tma_id` doesn't exist, that indexing expression
# evaluates to `character(0)`, so `matrix[character(0), ]` produces a
# 0-row matrix -- this script currently cannot produce a non-empty heatmap
# even once its input data exists, matching legacy's own literal behavior
# exactly.
#
# Reads $EXPORT_DIR/clinical.parquet, LEGACY_DATA_DIR's
# niche_frequencies_per_tma_id.parquet (currently missing, see above), and
# resources/colormaps.yaml. Writes to $OUTPUT_FIGURES_DIR/figure6/.

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
    # tma_id
  ) %>%
  distinct()

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
    annotation_colors[[col]] <- custom_annotation_colors[[col]]
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
ha <- rowAnnotation(df = df_metadata, col = annotation_colors, show_legend = c(FALSE, rep(TRUE, length(df_metadata) - 2), FALSE))

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

# reorder matrix columns (optional, consistent with YAML order)
matrix <- matrix[, intersect(names(niche_colors), colnames(matrix))]

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
