# Reproduce Supplementary Figure 5b: Kaplan-Meier overall survival by
# inflammation status (patient-level clinical variable, presence of at
# least one inflamed core).
#
# 1:1 port of 000_paper/11_niches/113_survival/risk_groups_label.R -- the
# "INFLAMMATION" section only (patient-level clinical `inflammation`
# variable, not niche-based -- distinct from Figure 6e's niche-derived risk
# score, figure6e_immune_risk_score_km.R). Not the patient-cluster (P1-P6),
# stromogenic (figureS6a_stromogenic_km.R), or Gleason-concordance
# (already figureS1b/figureS3b) sections of that same multi-analysis
# legacy script.
#
# Patient-level aggregation: max() of the per-core binary inflammation flag
# -- a patient counts as "inflamed" if at least one of their cores does.
#
# Disclosed fix: both ggsave() calls are commented out in legacy
# (computed-but-never-saved) -- enabled here using the exact args legacy's
# own commented-out calls specify (`plot = p$plot`, same width/height).
#
# Reads $EXPORT_DIR/clinical.parquet. Writes to
# $OUTPUT_FIGURES_DIR/figureS5/.

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

save_dir <- file.path(output_figures_dir, "figureS5")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))

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
plot_path <- file.path(save_dir, "figureS5b_progr_inflammation.pdf")
ggsave(plot_path, p_prog_inflam$plot, width = 6, height = 4)

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
plot_path <- file.path(save_dir, "figureS5b_survival_inflammation.pdf")
ggsave(plot_path, p_survival_inflam$plot, width = 6, height = 4)

cat("Saved Supplementary Figure 5b panels to", save_dir, "\n")
