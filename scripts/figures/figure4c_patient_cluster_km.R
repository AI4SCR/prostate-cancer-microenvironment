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

save_dir <- file.path(output_figures_dir, "figure4")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()

group.path <- file.path(legacy_dir, "5-niches", "barplot_data", "metadata_with_dendrogram_colors_label_pat_id.parquet")
df_patient <- read_parquet(group.path)

df <- df_patient |> filter(leaf_color_group != "black")

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
strata_order <- names(fit$strata)
strata_order <- gsub("cluster_group=", "", strata_order)
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

df$overall_survival <- ifelse(df$os_status == "alive", 0, 1)
fit <- survfit(Surv(last_fu, overall_survival) ~ cluster_group, data = df)

# Extract order used internally by survfit
strata_order <- names(fit$strata)
strata_order <- gsub("cluster_group=", "", strata_order)
names(custom_palette) <- names(fit$strata)
p_survival <- ggsurvplot(
  fit,
  data = df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = custom_palette,
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Group",
  risk.table.height = 0.25
)
p_survival_combined <- p_survival$plot / p_survival$table
plot_path <- file.path(save_dir, "figure4c_survival_by_patient_cluster.png")
ggsave(plot_path, p_survival_combined, width = 10, height = 6)

cat("Saved Figure 4c panels to", save_dir, "\n")
