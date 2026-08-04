# Reproduce Figure 4c: Kaplan-Meier survival by patient cluster.
#
# 1:1 port of 000_paper/sync_paper/03_survival/patient_risk_group_km.R
# (newly-pulled). Produces TWO plots, matching legacy exactly: progression-
# free survival (Surv(disease_progr_time, disease_progr)) and overall
# survival (Surv(last_fu, os_status=="dead")), both stratified by patient
# cluster (leaf_color_group), excluding "black" (dendrogram leaves above
# the color threshold, not a real cluster).
#
# Previously this was incorrectly combined into figure4_survival.R alongside
# the unrelated Cox PH panels (4d-e) -- that script's actual legacy source
# (survival-proportions.r) has no patient-cluster KM logic at all, only an
# ad hoc single-marker threshold-sweep KM exploration unrelated to any
# published panel. Split out into its own faithful port.
#
# Key correction vs. the earlier combined version: legacy uses an EXPLICIT
# custom color palette from the precomputed file's own `leaf_color` column
# (the actual dendrogram leaf colors, matching Figure 4a), not a default
# color scheme, and disables the confidence-interval band (`conf.int =
# FALSE`) -- the earlier version used default colors and added a CI band.
#
# survminer::ggsurvplot -> ggsurvfit substitution (already established
# elsewhere in this repo; survminer fails to compile in this environment,
# see REPRODUCIBILITY.md): ggsurvfit() + scale_color_manual() for the custom
# palette, add_risktable() for the risk table, no add_confidence_interval()
# call (matching conf.int = FALSE).
#
# Reads scripts/data/export.py's clinical.parquet and LEGACY_DATA_DIR's
# precomputed metadata_with_dendrogram_colors_label_pat_id.parquet (already
# carries os_status/disease_progr; only disease_progr_time/last_fu are
# joined in from clinical.parquet, matching legacy's left_join). Writes to
# OUTPUT_FIGURES_DIR/figure4/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(survival)
library(ggsurvfit)
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

df_patient <- read_parquet(file.path(
  legacy_dir, "5-niches", "barplot_data", "metadata_with_dendrogram_colors_label_pat_id.parquet"
))

df <- df_patient |> filter(leaf_color_group != "black")

df_colors <- df |>
  select(leaf_color_group, leaf_color) |>
  distinct() |>
  arrange(leaf_color_group)
custom_palette <- setNames(df_colors$leaf_color, df_colors$leaf_color_group)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
clinical_time <- clinical |>
  select(pat_id, disease_progr_time, last_fu) |>
  distinct()

df <- df |> left_join(clinical_time, by = "pat_id")
df$cluster_group <- factor(df$leaf_color_group, levels = names(custom_palette))

# %% Progression-free survival
# Legacy computes and displays this plot but never saves it -- its ggsave
# call is commented out, and references an undefined `result_dir` variable
# (a bug/incomplete cleanup left in the source). Matched verbatim: computed,
# not written to disk. This suggests the published Figure 4c is the overall
# survival panel below, not this one (consistent with REPRODUCIBILITY.md's
# note that progression-free appears in Supplementary Fig 3c instead).
fit_prog <- survfit2(Surv(disease_progr_time, disease_progr) ~ cluster_group, data = df)
p_prog <- fit_prog |>
  ggsurvfit() +
  scale_color_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  labs(title = "Figure 4c -- progression-free survival by patient cluster", x = "Time", y = "Progression-free survival probability") +
  add_risktable()
p_prog

# %% Overall survival
df$overall_survival <- ifelse(df$os_status == "alive", 0, 1)
fit_survival <- survfit2(Surv(last_fu, overall_survival) ~ cluster_group, data = df)
p_survival <- fit_survival |>
  ggsurvfit() +
  scale_color_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  labs(title = "Figure 4c -- overall survival by patient cluster", x = "Time", y = "Survival probability") +
  add_risktable()
ggsave(file.path(save_dir, "figure4c_survival_by_patient_cluster.png"), p_survival, width = 10, height = 6, dpi = 200)

cat("Saved Figure 4c panels to", save_dir, "\n")
