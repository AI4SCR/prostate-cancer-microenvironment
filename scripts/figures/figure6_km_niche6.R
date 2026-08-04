# Reproduce Figure 6b: KM survival/progression by high-vs-low niche
# abundance (binary median split), scanned across every niche -- niche 6 is
# the one that makes the figure, but the old script (and this port) produces
# the same panel for all niches, letting you pick niche 6 out of the output.
#
# 1:1 port of the old repo's
# 000_paper/11_niches/113_survival/niche_kaplan_meier_binary.R (paths only
# changed, unused `ggpubr` import dropped; see figure_script_mapping.md).
# clusters_annotated_v2.parquet has no reproducing script in this repo
# (figure5_niche_annotation.py is blocked -- see that script's docstring),
# so it's read from LEGACY_DATA_DIR's precomputed copy; clinical.parquet
# comes from EXPORT_DIR.
#
# Plotting uses `ggsurvfit`/`survfit2` instead of the original's
# `survminer::ggsurvplot` -- survminer's dependency chain (via ggpubr ->
# rstatix -> car -> pbkrtest -> doBy -> Deriv) fails to compile against this
# R 4.4.1 install (Deriv's C++ source uses `R_ClosureFormals`, not available
# here). `ggsurvfit` is already the established KM-plotting convention in
# this repo (see figure4_survival.R) and is installed; the statistical
# analysis (survfit/survdiff/coxph) is unchanged.
#
# Writes one KM PDF per niche (survival + progression) to
# $EXPORT_DIR/figures/figure6/niche_km/threshold_50/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(rlang)

export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

save_dir <- file.path(output_figures_dir, "figure6", "niche_km")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()

compute_label_frequency <- function(data, level, pseudocount = 1) {
  level_sym <- rlang::sym(level)
  data_summary <- data %>%
    dplyr::group_by(sample_name, !!level_sym) %>%
    dplyr::summarise(count = dplyr::n() + pseudocount, .groups = "drop") %>%
    tidyr::complete(sample_name, !!level_sym, fill = list(count = pseudocount))
  data_summary %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(proportion = count / sum(count)) %>%
    dplyr::ungroup()
}

df_clusters <- read_parquet(file.path(legacy_dir, "5-niches", "annotation", "clusters_annotated_v2.parquet"))
df_clusters[["sample_name"]] <- df_clusters[["tma_id"]]

df_freqs <- compute_label_frequency(df_clusters, level = "niche", pseudocount = 0)
df_props <- df_freqs %>%
  select(tma_id = sample_name, niche, proportion) %>%
  pivot_wider(names_from = niche, values_from = proportion, values_fill = 0)

cols <- colnames(df_props)[2:ncol(df_props)]
qs <- c(0.25, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)

quantiles <- list()
for (col in cols) {
  x_nz <- df_props[[col]][df_props[[col]] > 0]
  quantiles[[col]] <- quantile(x_nz, probs = qs, na.rm = TRUE)
}

threshold <- "50%"
thr <- sapply(cols, function(col) quantiles[[col]][[threshold]])
names(thr) <- cols

save_dir <- file.path(save_dir, paste0("threshold_", substr(threshold, 1, 2)))
progression_dir <- file.path(save_dir, "progression_free_survival")
survival_dir <- file.path(save_dir, "overall_survival")
dir.create(progression_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(survival_dir, recursive = TRUE, showWarnings = FALSE)

df_binary <- df_props %>%
  mutate(across(-tma_id, ~ ifelse(.x >= thr[cur_column()], 1, 0)))

clinical$os_status <- ifelse(clinical$os_status == "dead", 1, 0)
progression <- clinical %>% select(pat_id, disease_progr, disease_progr_time) %>% distinct()
death <- clinical %>% select(pat_id, os_status, last_fu) %>% distinct()

results_os <- list()
results_prog <- list()

for (col in cols) {
  df_niche <- df_binary
  df_niche[["target"]] <- df_niche[[col]]
  df_patient <- clinical %>%
    select(pat_id, tma_id) %>%
    distinct() %>%
    inner_join(df_niche, by = "tma_id") %>%
    select(pat_id, risk_group = target) %>%
    group_by(pat_id) %>%
    summarise(risk_group = max(risk_group), .groups = "drop")

  df_analysis <- df_patient %>%
    inner_join(progression, by = "pat_id") %>%
    inner_join(death, by = "pat_id")

  fit <- survfit2(Surv(last_fu, os_status) ~ risk_group, data = df_analysis)
  form <- as.formula(fit$call$formula)
  sd <- survdiff(form, data = df_analysis)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  cox <- coxph(formula = form, data = df_analysis)
  s_cox <- summary(cox)$coefficients[, c("exp(coef)", "Pr(>|z|)")]
  results_os[[col]] <- list(pval = pval, cox_coef = s_cox[1], cox_pval = s_cox[2])
  p1 <- fit |>
    ggsurvfit() +
    labs(title = paste("Survival by", col, "high vs low"), x = "Time", y = "Survival probability") +
    add_risktable() +
    add_pvalue()
  pdf(file.path(survival_dir, paste0("km_survival_os_status_", col, ".pdf")), width = 8, height = 6)
  print(p1)
  dev.off()

  fit_prog <- survfit2(Surv(disease_progr_time, disease_progr) ~ risk_group, data = df_analysis)
  form <- as.formula(fit_prog$call$formula)
  sd <- survdiff(form, data = df_analysis)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  cox <- coxph(formula = form, data = df_analysis)
  s_cox <- summary(cox)$coefficients[, c("exp(coef)", "Pr(>|z|)")]
  results_prog[[col]] <- list(pval = pval, cox_coef = s_cox[1], cox_pval = s_cox[2])
  p2 <- fit_prog |>
    ggsurvfit() +
    labs(title = paste("Progression-free by", col, "high vs low"), x = "Time", y = "Progression-free probability") +
    add_risktable() +
    add_pvalue()
  pdf(file.path(progression_dir, paste0("km_progression_disease_progr_", col, ".pdf")), width = 8, height = 6)
  print(p2)
  dev.off()
}

df_os <- do.call(rbind, lapply(names(results_os), function(col) {
  data.frame(niche = col, pval = results_os[[col]]$pval, cox_coef = results_os[[col]]$cox_coef, cox_pval = results_os[[col]]$cox_pval)
}))
df_os$qval <- round(p.adjust(df_os$pval, method = "BH"), 4)
df_os$cox_qval <- round(p.adjust(df_os$cox_pval, method = "BH"), 4)
write.csv(df_os, file.path(save_dir, "figure6b_niche_km_os_summary.csv"), row.names = FALSE)

df_prog <- do.call(rbind, lapply(names(results_prog), function(col) {
  data.frame(niche = col, pval = results_prog[[col]]$pval, cox_coef = results_prog[[col]]$cox_coef, cox_pval = results_prog[[col]]$cox_pval)
}))
df_prog$qval <- round(p.adjust(df_prog$pval, method = "BH"), 4)
df_prog$cox_qval <- round(p.adjust(df_prog$cox_pval, method = "BH"), 4)
write.csv(df_prog, file.path(save_dir, "figure6b_niche_km_progression_summary.csv"), row.names = FALSE)

cat("Saved per-niche KM panels + summary tables to", save_dir, "\n")
