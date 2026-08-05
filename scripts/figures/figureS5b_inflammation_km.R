library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(survival)
library(survminer)
library(patchwork)

export_dir <- Sys.getenv("EXPORT_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))

figures_dir <- file.path(output_figures_dir, "figureS5")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))

dir_inflam <- file.path(figures_dir, "inflammation")
df_inflam <- clinical %>%
  select(pat_id, inflammation) %>%
  filter(!is.na(inflammation)) %>%
  mutate(
    inflammation_bin = case_when(
      inflammation == "yes" ~ 1,
      inflammation == "no" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(pat_id) %>%
  summarise(
    inflammation = max(inflammation_bin, na.rm = TRUE)
  ) %>%
  ungroup()

df_outcome <- clinical %>%
  select(pat_id, disease_progr, disease_progr_time, last_fu, os_status) %>%
  distinct()

df_final <- df_inflam %>%
  left_join(df_outcome, by = "pat_id")
df_final$overall_survival <- ifelse(df_final$os_status == "alive", 0, 1)

fit <- survfit(Surv(disease_progr_time, disease_progr) ~ inflammation, data = df_final)
p_prog_inflam <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "Inflammation",
  risk.table.height = 0.25
)
plot_name <- "progr_inflammation.pdf"
plot_path <- file.path(dir_inflam, plot_name)
ggsave(plot_path, p_prog_inflam$plot, width = 6, height = 4)
p_prog_inflam$plot
fit <- survfit(Surv(last_fu, overall_survival) ~ inflammation, data = df_final)
p_survival_inflam <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Inflammation",
  risk.table.height = 0.25
)
p <- p_survival_inflam$plot / p_survival_inflam$table
plot_name <- "survival_inflammation_with_table.pdf"
plot_path <- file.path(dir_inflam, plot_name)
ggsave(plot_path, p, width = 6, height = 4)
p_survival_inflam$plot
