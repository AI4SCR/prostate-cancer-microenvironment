library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

export_dir <- Sys.getenv("EXPORT_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))

figures_dir <- file.path(output_figures_dir, "figureS1")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))

filter_ <- duplicated(clinical$tma_id)
tmas <- clinical[!filter_, ]

var_name <- "gleason_grp" # x-axis categories in the heatmap
order_name <- "gs_grp" # patient-level ordering variable

# 1) patient order (unique gs_grp per patient)
pat_order <- tmas %>%
  select(pat_id, !!sym(order_name)) %>%
  filter(!is.na(.data[[order_name]]), .data[[order_name]] != "nan") %>%
  distinct() %>%
  arrange(as.numeric(.data[[order_name]]), pat_id) %>%
  mutate(pat_id = factor(pat_id, levels = pat_id))

# 2) counts long table
pdat <- as.data.frame.matrix(table(tmas[[var_name]], tmas$pat_id)) |>
  rownames_to_column(var = var_name) |>
  pivot_longer(
    cols = -all_of(var_name),
    names_to = "pat_id",
    values_to = "value"
  ) |>
  semi_join(pat_order, by = "pat_id") |>
  mutate(pat_id = factor(pat_id, levels = levels(pat_order$pat_id)))

# 3) main heatmap
g_main <- ggplot(
  pdat,
  aes(
    y = pat_id,
    x = factor(.data[[var_name]]),
    fill = factor(value, levels = 0:4)
  )
) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_manual(
    values = c("0" = "white", "4" = "#fde725", "3" = "#5ec962", "2" = "#21918c", "1" = "#440154"),
    drop = FALSE,
    name = "Count"
  ) +
  labs(y = "Patient", x = var_name) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )

gs_palette <- c(
  "1" = "#f6d2d2",
  "2" = "#f1abab",
  "3" = "#f96b6b",
  "4" = "#f44336",
  "5" = "#b71c1c",
  "nan" = "#9e9e9e"
)

g_strip <- pat_order %>%
  mutate(
    gs_grp = as.character(.data[[order_name]]),
    gs_grp = factor(gs_grp, levels = names(gs_palette)),
    x = "gs_grp"
  ) %>%
  ggplot(aes(x = x, y = pat_id, fill = gs_grp)) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_manual(values = gs_palette, drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "GS group") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

p <- (g_strip + g_main) + plot_layout(widths = c(1, 10))
plot_name <- "heatmap_gleason_grp_by_gs_grp.pdf"
plot_path <- file.path(figures_dir, plot_name)

ggsave(plot_path, p, width = 10, height = 20)
cat("Saved Supplementary Figure 1b to", plot_path, "\n")
