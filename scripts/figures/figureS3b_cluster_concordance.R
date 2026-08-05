library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

figures_dir <- file.path(output_figures_dir, "figureS3")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()
cat("Number of patients:", num.patients, "\n")

## read patient-level and core-level cluster group
barplot_data_dir <- file.path(legacy_dir, "5-niches", "barplot_data")

path_patient <- file.path(barplot_data_dir, "metadata_with_dendrogram_colors_label_pat_id.parquet")
df_patient <- read_parquet(path_patient)
path_tma <- file.path(barplot_data_dir, "metadata_with_dendrogram_colors_label_tma_id.parquet")
df_tma <- read_parquet(path_tma)

df_patient <- df_patient %>%
  select(pat_id, leaf_color_group, leaf_color)
df_tma <- df_tma %>%
  select(tma_id, leaf_color_group, leaf_color, pat_id)
df_tma <- df_tma %>%
  rename(
    pat_id = pat_id,
    cluster_group_tma = leaf_color_group,
    leaf_color_tma = leaf_color
  )
df_patient <- df_patient %>%
  rename(
    cluster_group_patient = leaf_color_group,
    leaf_color_patient = leaf_color
  )

df_tma <- df_tma %>%
  left_join(df_patient, by = "pat_id")

var_name <- "cluster_group_tma" # x-axis categories in the heatmap
order_name <- "cluster_group_patient" # patient-level ordering variable

df_tma[[order_name]] <- as.factor(df_tma[[order_name]])
# 1) patient order (unique gs_grp per patient)
pat_order <- df_tma %>%
  select(pat_id, !!sym(order_name)) %>%
  filter(!is.na(.data[[order_name]]), .data[[order_name]] != "nan") %>%
  distinct() %>%
  arrange(as.numeric(.data[[order_name]]), pat_id) %>%
  mutate(pat_id = factor(pat_id, levels = pat_id))

# 2) counts long table
pdat <- as.data.frame.matrix(table(df_tma[[var_name]], df_tma$pat_id)) |>
  rownames_to_column(var = var_name) |>
  pivot_longer(
    cols = -all_of(var_name),
    names_to = "pat_id",
    values_to = "value"
  ) |>
  semi_join(pat_order, by = "pat_id") |>
  mutate(pat_id = factor(pat_id, levels = levels(pat_order$pat_id)))

tma_palette <- df_tma %>%
  select(cluster_group_tma, leaf_color_tma) %>%
  distinct() %>%
  arrange(cluster_group_tma)

pat_palette <- df_patient %>%
  select(cluster_group_patient, leaf_color_patient) %>%
  distinct() %>%
  arrange(cluster_group_patient)

# 3) main heatmap
# palette vectors
tma_cols <- setNames(
  tma_palette$leaf_color_tma,
  tma_palette$cluster_group_tma
)

g_main <- ggplot(
  pdat,
  aes(
    y = pat_id,
    x = factor(.data[[var_name]]),
    fill = factor(value, levels = 0:4)
  )
) +
  geom_tile(color = "black", linewidth = 0.1) +

  # dummy points only for legend
  geom_point(
    data = data.frame(cluster_group_tma = names(tma_cols)),
    aes(x = NA, y = NA, color = cluster_group_tma),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("0" = "white", "4" = "#fde725", "3" = "#5ec962", "2" = "#21918c", "1" = "#440154"),
    drop = FALSE,
    name = "Count"
  ) +
  scale_color_manual(
    values = tma_cols,
    name = "TMA cluster"
  ) +
  labs(y = "Patient", x = var_name) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )

## make to list
pat_palette <- setNames(pat_palette$leaf_color_patient, pat_palette$cluster_group_patient)

g_strip <- pat_order %>%
  mutate(
    gs_grp = as.character(.data[[order_name]]),
    gs_grp = factor(gs_grp, levels = names(pat_palette)),
    x = "gs_grp"
  ) %>%
  ggplot(aes(x = x, y = pat_id, fill = gs_grp)) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_manual(values = pat_palette, drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "GS group") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

(g_strip + g_main) + plot_layout(widths = c(1, 10))

# save as pdf
plot_name <- "heatmap_cluster_group_tma_by_patient_cluster_group.pdf"
plot_path <- file.path(save_dir, plot_name)
ggsave(plot_path, (g_strip + g_main) + plot_layout(widths = c(1, 10)), width = 10, height = 20)
cat("Saved Supplementary Figure 3b to", plot_path, "\n")
