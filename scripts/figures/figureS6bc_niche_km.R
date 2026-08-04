# Reproduce Supplementary Figure 6b-c: KM progression-free survival by
# high-vs-low niche 2 / niche 8 abundance (binary median split).
#
# 1:1 port of the old repo's
# 000_paper/11_niches/113_survival/niche_kaplan_meier_binary.R, trimmed to
# niches 2 and 8 only -- the full 18-niche loop is already the
# already-validated figure6_km_niche6.R (not touched here). Niche 2 =
# `luminal_infiltrated`, niche 8 = `tumor_CAF1_lymphocytes`, per the fixed
# niche_order list used throughout the legacy niche-visualization scripts
# (see figureS4b_niche_mean_composition.py's NICHE_ORDER) -- this mapping is
# not stored in any data file and is only independently confirmed for
# niches 6 and 16-18; flagged in open-questions.md. Bundled into one script
# since both panels are the same computation for two niches within one
# supplementary figure (same precedent as Figure 2b's compartment-subset
# loop, not a "different code path per panel").
#
# clusters_annotated_v2.parquet has no reproducing script in this repo, so
# it's read from LEGACY_DATA_DIR's precomputed copy; clinical.parquet comes
# from EXPORT_DIR.
#
# Writes KM PDFs (survival + progression) to
# $OUTPUT_FIGURES_DIR/figureS6/niche_km/threshold_50/.

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

save_dir <- file.path(output_figures_dir, "figureS6", "niche_km")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

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

df_freqs <- compute_label_frequency(df_clusters, level = "niche", pseudocount = 0)
df_props <- df_freqs %>%
  select(tma_id = sample_name, niche, proportion) %>%
  pivot_wider(names_from = niche, values_from = proportion, values_fill = 0)

cols <- c("luminal_infiltrated", "tumor_CAF1_lymphocytes") # niches 2 and 8
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
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

df_binary <- df_props %>%
  mutate(across(all_of(cols), ~ ifelse(.x >= thr[cur_column()], 1, 0)))

clinical$os_status <- ifelse(clinical$os_status == "dead", 1, 0)
progression <- clinical %>% select(pat_id, disease_progr, disease_progr_time) %>% distinct()
death <- clinical %>% select(pat_id, os_status, last_fu) %>% distinct()

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

  fit <- survfit(Surv(last_fu, os_status) ~ risk_group, data = df_analysis)
  p1 <- ggsurvplot(
    fit,
    data = df_analysis,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    palette = "Set2",
    xlab = "Time",
    ylab = "Survival probability",
    legend.title = "Group",
    risk.table.height = 0.25,
    title = paste("Survival by", col, "high vs low")
  )
  pdf(file.path(save_dir, paste0("km_survival_os_status_", col, ".pdf")), width = 8, height = 6)
  print(p1)
  dev.off()

  fit_prog <- survfit(Surv(disease_progr_time, disease_progr) ~ risk_group, data = df_analysis)
  p2 <- ggsurvplot(
    fit_prog,
    data = df_analysis,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    palette = "Set2",
    xlab = "Time",
    ylab = "Progression-free probability",
    legend.title = "Group",
    risk.table.height = 0.25,
    title = paste("Progression-free by", col, "high vs low")
  )
  pdf(file.path(save_dir, paste0("km_progression_disease_progr_", col, ".pdf")), width = 8, height = 6)
  print(p2)
  dev.off()
}

cat("Saved Supplementary Figure 6b-c (niches 2, 8) KM panels to", save_dir, "\n")
