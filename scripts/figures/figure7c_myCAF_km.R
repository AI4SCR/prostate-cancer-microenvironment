# Reproduce Figure 7c: KM progression-free survival by high-vs-low myCAF
# (CD105-high CAF1) abundance (binary 75th-percentile split).
#
# 1:1 port of the old repo's
# 000_paper/11_niches/113_survival/label_kaplan_meier_binary.R, trimmed to
# the single `stromal-CAF1(CD105+)` label only -- remove the loop over all
# ~34 other cell-type labels, they belong to no current panel. "myCAF" =
# `stromal-CAF1(CD105+)`: the paper text states "CD105high are annotated as
# myCAFs", and `resources/colormaps.yaml`'s `label:` key confirms
# `stromal-CAF1(CD105+)` is the CD105-high CAF1 variant (not
# `stromal-CAF2(AR+)`, an earlier incorrect guess -- see open-questions.md).
#
# Same structure as niche_kaplan_meier_binary.R (already ported as
# figure6_km_niche6.R/figure7b_niche9_km.R/figureS6bc_niche_km.R), but with
# level="label" instead of "niche", threshold=75th percentile (not 50th),
# and `conf.int = TRUE` in both ggsurvplot() calls (legacy's own, not a
# deviation -- differs from the niche-level scripts' `conf.int = FALSE`).
#
# Disclosed fix: both ggsave() calls are commented out in legacy
# (computed-but-never-saved) -- enabled here using `p$plot` (matching
# legacy's own commented-out call args).
#
# clusters_annotated_v2.parquet has no reproducing script in this repo, so
# it's read from LEGACY_DATA_DIR's precomputed copy; clinical.parquet comes
# from EXPORT_DIR.
#
# Writes KM PDFs (survival + progression) to
# $OUTPUT_FIGURES_DIR/figure7/label_km/threshold_75/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(rlang)

export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

save_dir <- file.path(output_figures_dir, "figure7", "label_km")

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))

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

df_freqs <- compute_label_frequency(df_clusters, level = "label", pseudocount = 0)
df_props <- df_freqs %>%
  select(tma_id = sample_name, label, proportion) %>%
  pivot_wider(names_from = label, values_from = proportion, values_fill = 0)

cols <- "stromal-CAF1(CD105+)" # myCAF
qs <- c(0.25, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)

quantiles <- list()
for (col in cols) {
  x_nz <- df_props[[col]][df_props[[col]] > 0]
  quantiles[[col]] <- quantile(x_nz, probs = qs, na.rm = TRUE)
}

threshold <- "75%"
thr <- sapply(cols, function(col) quantiles[[col]][[threshold]])
names(thr) <- cols

save_dir <- file.path(save_dir, paste0("threshold_", substr(threshold, 1, 2)))
progression_dir <- file.path(save_dir, "progression_free_survival")
survival_dir <- file.path(save_dir, "overall_survival")
dir.create(progression_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(survival_dir, recursive = TRUE, showWarnings = FALSE)

df_binary <- df_props %>%
  mutate(across(all_of(cols), ~ ifelse(.x >= thr[cur_column()], 1, 0)))

clinical$os_status <- ifelse(clinical$os_status == "dead", 1, 0)
progression <- clinical %>% select(pat_id, disease_progr, disease_progr_time) %>% distinct()
death <- clinical %>% select(pat_id, os_status, last_fu) %>% distinct()

for (col in cols) {
  df_label <- df_binary
  df_label[["target"]] <- df_label[[col]]
  df_patient <- clinical %>%
    select(pat_id, tma_id) %>%
    distinct() %>%
    inner_join(df_label, by = "tma_id") %>%
    select(pat_id, risk_group = target) %>%
    group_by(pat_id) %>%
    summarise(risk_group = max(risk_group), .groups = "drop")

  df_analysis <- df_patient %>%
    inner_join(progression, by = "pat_id") %>%
    inner_join(death, by = "pat_id")

  fit <- survfit(Surv(last_fu, os_status) ~ risk_group, data = df_analysis)
  p1 <- ggsurvplot(
    fit,
    data = df_analysis,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    palette = "Set2",
    xlab = "Time",
    ylab = "Survival probability",
    legend.title = "Group",
    risk.table.height = 0.25,
    title = paste("Survival by", col, "high vs low")
  )
  plot_path <- file.path(survival_dir, paste0("km_survival_os_status_", col, ".pdf"))
  ggsave(plot_path, p1$plot, width = 8, height = 6, dpi = 300)

  fit_prog <- survfit(Surv(disease_progr_time, disease_progr) ~ risk_group, data = df_analysis)
  p2 <- ggsurvplot(
    fit_prog,
    data = df_analysis,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    palette = "Set2",
    xlab = "Time",
    ylab = "Progression-free probability",
    legend.title = "Group",
    risk.table.height = 0.25,
    title = paste("Progression-free by", col, "high vs low")
  )
  plot_path <- file.path(progression_dir, paste0("km_progression_disease_progr_", col, ".pdf"))
  ggsave(plot_path, p2$plot, width = 8, height = 6, dpi = 300)
}

cat("Saved Figure 7c (myCAF) KM panels to", save_dir, "\n")
