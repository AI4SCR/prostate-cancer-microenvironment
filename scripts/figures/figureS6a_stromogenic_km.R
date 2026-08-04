# Reproduce Supplementary Figure 6a: Kaplan-Meier progression-free survival
# by stromogenic status (patient-level clinical variable, presence of at
# least one histologically stromogenic core).
#
# 1:1 port of 000_paper/11_niches/113_survival/risk_groups_label.R -- the
# "STROMOGENIC" section only (`stromogenic_smc_loss_reactive_stroma_present`).
# Not the patient-cluster, inflammation (figureS5b_inflammation_km.R), or
# Gleason-concordance sections of that same multi-analysis legacy script.
#
# Patient-level aggregation: max() of the per-core binary stromogenic flag
# -- same pattern as figureS5b_inflammation_km.R's inflammation section.
#
# Disclosed fix: both ggsave() calls are commented out in legacy
# (computed-but-never-saved) -- enabled here using the exact args legacy's
# own commented-out calls specify (`plot = p$plot`, same width/height).
#
# Reads $EXPORT_DIR/clinical.parquet. Writes to
# $OUTPUT_FIGURES_DIR/figureS6/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(survival)
library(survminer)

export_dir <- Sys.getenv("EXPORT_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))

save_dir <- file.path(output_figures_dir, "figureS6")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))

df_stromo <- clinical %>%
  select(pat_id, stromogenic_smc_loss_reactive_stroma_present) %>%
  filter(!is.na(stromogenic_smc_loss_reactive_stroma_present)) %>%
  mutate(
    stromogenic_bin = case_when(
      stromogenic_smc_loss_reactive_stroma_present == "yes" ~ 1,
      stromogenic_smc_loss_reactive_stroma_present == "no" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(pat_id) %>%
  summarise(
    stromogenic = max(stromogenic_bin, na.rm = TRUE)
  ) %>%
  ungroup()

df_outcome <- clinical %>%
  select(pat_id, disease_progr, disease_progr_time, last_fu, os_status) %>%
  distinct()

df_final <- df_stromo %>%
  left_join(df_outcome, by = "pat_id")
df_final$overall_survival <- ifelse(df_final$os_status == "alive", 0, 1)

fit <- survfit(Surv(disease_progr_time, disease_progr) ~ stromogenic, data = df_final)
p_prog_stromo <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "stromogenic",
  risk.table.height = 0.25
)
plot_path <- file.path(save_dir, "figureS6a_progr_stromogenic.pdf")
ggsave(plot_path, p_prog_stromo$plot, width = 6, height = 4)

fit <- survfit(Surv(last_fu, overall_survival) ~ stromogenic, data = df_final)
p_survival_stromo <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "stromogenic",
  risk.table.height = 0.25
)
plot_path <- file.path(save_dir, "figureS6a_survival_stromogenic.pdf")
ggsave(plot_path, p_survival_stromo$plot, width = 6, height = 4)

cat("Saved Supplementary Figure 6a panels to", save_dir, "\n")
