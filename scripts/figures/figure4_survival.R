# Reproduce Figure 4c-e: survival analysis on patient-level cell-type composition.
#
# 4c: Kaplan-Meier overall survival stratified by patient cluster, two-sided
#     log-rank test.
# 4d-e: univariate Cox PH regression per cell type (CLR-transformed
#     max-pooled patient-level proportions) for overall survival and disease
#     progression, BH-adjusted p-values, plotted as a hazard-ratio forest plot.
#
# Figure 4c's cluster labels come from LEGACY_DATA_DIR's precomputed
# metadata_with_dendrogram_colors_label_pat_id.parquet, NOT from
# figure4_patient_clustering.py's own output. That script's reconstruction
# (JSD linkage + scipy fcluster(criterion="distance")) is now confirmed
# correct in its core computation -- its 4 largest clusters (93/63/11/9
# patients) are byte-identical in size to this precomputed file -- but
# fcluster's flat-clustering heuristic assigns the smallest, near-singleton
# leaves differently than matplotlib's dendrogram(color_threshold=...)
# leaf-coloring (the actual algorithm that produced this file) does. Rather
# than accept that residual mismatch, this reads the real precomputed
# partition directly. See REPRODUCIBILITY.md's Figure 4a/c note.
#
# Reads scripts/00-data-export/export_for_r.py's clinical.parquet,
# figure4_patient_clustering.py's composition output (for 4d-e, unaffected
# by the clustering question), and the precomputed cluster file above.
# Writes to EXPORT_DIR/figures/figure4/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(survival)
library(ggsurvfit)
library(compositions)
library(ggplot2)

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

composition <- read_parquet(file.path(save_dir, "figure4a_patient_composition.parquet"))
patient_clusters <- read_parquet(file.path(
  legacy_dir, "5-niches", "barplot_data", "metadata_with_dendrogram_colors_label_pat_id.parquet"
)) |>
  transmute(pat_id, patient_cluster = leaf_color_group)
clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))

patient_clinical <- clinical |>
  select(pat_id, last_fu, os_status, disease_progr, disease_progr_time) |>
  distinct(pat_id, .keep_all = TRUE)
stopifnot("expected one clinical row per patient" = !any(duplicated(patient_clinical$pat_id)))

# %% Figure 4c: KM by patient cluster
km_data <- patient_clusters |>
  inner_join(patient_clinical, by = "pat_id") |>
  mutate(event = os_status == "dead", time = last_fu)
stopifnot("missing survival time/event for some patients" = !any(is.na(km_data$time) | is.na(km_data$event)))

fit <- survfit2(Surv(time, event) ~ patient_cluster, data = km_data)
logrank <- survdiff(Surv(time, event) ~ patient_cluster, data = km_data)
p_value <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)

p <- fit |>
  ggsurvfit() +
  labs(
    title = sprintf("Figure 4c -- overall survival by patient cluster (log-rank p=%.3g)", p_value),
    x = "Months", y = "Overall survival probability"
  ) +
  add_confidence_interval() +
  add_risktable()
ggsave(file.path(save_dir, "figure4c_km_patient_clusters.png"), p, width = 9, height = 8, dpi = 200)

# %% Figure 4d-e: Cox PH per cell type, CLR-transformed proportions
cell_types <- setdiff(colnames(composition), "pat_id")
clr_mat <- clr(as.matrix(composition[, cell_types]))  # clr() needs the full composition row at once
clr_composition <- composition |>
  select(pat_id) |>
  bind_cols(as.data.frame(clr_mat))

cox_data <- clr_composition |> inner_join(patient_clinical, by = "pat_id")

fit_cox <- function(data, cell_type, event_col, time_col) {
  d <- data |> transmute(x = .data[[cell_type]], time = .data[[time_col]], event = .data[[event_col]] == ifelse(event_col == "os_status", "dead", 1))
  fit <- coxph(Surv(time, event) ~ x, data = d)
  s <- summary(fit)
  tibble(
    cell_type = cell_type,
    hr = unname(s$coefficients[1, "exp(coef)"]),
    hr_lower = s$conf.int[1, "lower .95"],
    hr_upper = s$conf.int[1, "upper .95"],
    p_value = unname(s$coefficients[1, "Pr(>|z|)"]),
  )
}

run_cox_panel <- function(event_col, time_col, label) {
  results <- purrr::map_dfr(cell_types, ~ fit_cox(cox_data, .x, event_col, time_col)) |>
    mutate(p_adj = p.adjust(p_value, method = "fdr")) |>
    arrange(p_value)
  write_csv <- utils::write.csv
  write_csv(results, file.path(save_dir, sprintf("figure4_cox_%s.csv", label)), row.names = FALSE)

  p <- results |>
    mutate(cell_type = factor(cell_type, levels = rev(cell_type))) |>
    ggplot(aes(x = hr, y = cell_type, color = p_adj < 0.05)) +
    geom_point() +
    geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), height = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    labs(title = sprintf("Figure 4%s -- Cox PH hazard ratios (%s)", label, event_col), x = "Hazard ratio (log scale)", y = NULL) +
    theme(legend.position = "bottom")
  ggsave(file.path(save_dir, sprintf("figure4_cox_%s.png", label)), p, width = 8, height = 10, dpi = 200)
  results
}

os_results <- run_cox_panel("os_status", "last_fu", "d")
progression_results <- run_cox_panel("disease_progr", "disease_progr_time", "e")

cat("Saved figure 4c-e panels to", save_dir, "\n")
