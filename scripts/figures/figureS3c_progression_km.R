library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(survival)
library(survminer)
library(patchwork)

base_dir <- Sys.getenv("BASE_DIR")
export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(export_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

figures_dir <- file.path(output_figures_dir, "figureS3")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()

group.path <- file.path(legacy_dir, "5-niches", "barplot_data", "metadata_with_dendrogram_colors_label_pat_id.parquet")
df_patient <- read_parquet(group.path)

df <- df_patient |> filter(leaf_color_group != "black")
df$leaf_color_group <- sub("^C", "P", df$leaf_color_group) # "C1".."C6" -> "P1".."P6", see docstring

df_colors <- df |>
  select(leaf_color_group, leaf_color) |>
  distinct() |>
  arrange(leaf_color_group)

custom_palette <- df_colors$leaf_color
names(custom_palette) <- df_colors$leaf_color_group

## merge with disease_progr, last_fu
clinical_time <- clinical |>
  select(pat_id, disease_progr_time, last_fu) |>
  distinct()

df <- df |> left_join(clinical_time, by = "pat_id")

df$cluster_group <- factor(df$leaf_color_group, levels = names(custom_palette))
fit <- survfit(Surv(disease_progr_time, disease_progr) ~ cluster_group, data = df)

# Extract order used internally by survfit
names(custom_palette) <- names(fit$strata)
p_prog <- ggsurvplot(
  fit,
  data = df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = custom_palette,
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "Group",
  risk.table.height = 0.25
)

p_prog$plot

plot_name <- "progr_patient_cluster_group_final_all.pdf"
plot_path <- file.path(result_dir, plot_name)
# ggsave(plot_path, p_prog$plot, width = 10, height = 6)
