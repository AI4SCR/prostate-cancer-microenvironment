# Reproduce Supplementary Figure 3c: Kaplan-Meier progression-free survival
# by patient cluster (P1-P6).
#
# 1:1 port of 000_paper/sync_paper/03_survival/patient_risk_group_km.R --
# the progression-free branch only. The overall-survival branch of that same
# script is Figure 4c's already-validated script,
# figure4c_patient_cluster_km.R, which is not touched here; this is an
# independent, standalone port of the same data-loading/relabeling logic
# (duplicated rather than shared, per this repo's own anti-pattern rule
# against cross-script imports of reusable logic).
#
# Disclosed fix: legacy computes this exact plot but never saves it -- its
# ggsave call is commented out and references an undefined `result_dir`
# variable (an authoring-artifact bug, matched verbatim as dead code in
# figure4c_patient_cluster_km.R's port of the overall-survival branch, since
# that branch's own save works fine). Since this script's entire purpose
# *is* this panel, the save is enabled here -- same precedent as
# figure6_niche_abundance_heatmap.R's already-enabled commented-out
# pdf()/dev.off().
#
# Cluster labels: "C1".."C6" -> "P1".."P6", same disclosed correction as
# figure4c_patient_cluster_km.R (see that script's docstring for the full
# rationale).
#
# Reads scripts/data/export.py's clinical.parquet and LEGACY_DATA_DIR's
# precomputed metadata_with_dendrogram_colors_label_pat_id.parquet. Writes
# to OUTPUT_FIGURES_DIR/figureS3/.

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

save_dir <- file.path(output_figures_dir, "figureS3")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

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

p_prog_combined <- p_prog$plot / p_prog$table
plot_path <- file.path(save_dir, "figureS3c_progression_by_patient_cluster.png")
ggsave(plot_path, p_prog_combined, width = 10, height = 6)

cat("Saved Supplementary Figure 3c to", plot_path, "\n")
