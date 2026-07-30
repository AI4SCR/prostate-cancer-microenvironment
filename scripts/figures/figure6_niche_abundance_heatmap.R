# Reproduce Figure 6a: per-patient niche-proportion heatmap (+ the
# dendrogram-split risk-group KM analysis the same script also produces).
#
# 1:1 port of the old repo's
# 000_paper/11_niches/111_heatmaps/patient_heatmap.R (paths only changed;
# see figure_script_mapping.md). Niche composition/annotation tables have
# no reproducing script in this repo, so they come from LEGACY_DATA_DIR;
# clinical.parquet comes from EXPORT_DIR.
#
# Plotting uses `ggsurvfit`/`survfit2` instead of the original's
# `survminer::ggsurvplot` -- survminer fails to compile in this environment
# (see figure6_km_niche6.R's docstring for the exact error). `ggsurvfit` is
# already this repo's established KM-plotting convention (figure4_survival.R).
#
# Writes to $EXPORT_DIR/figures/figure6/.

library(dotenv)
load_dot_env()

library(arrow)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(entropy)
library(survival)
library(ggsurvfit)

export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

save_dir <- file.path(export_dir, "figures", "figure6")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

niches_dir <- file.path(legacy_dir, "5-niches")

df_props <- read_parquet(file.path(niches_dir, "frequencies", "stacked_barplots", "props_niche_pat_id.parquet"))
clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()

matrix <- df_props %>%
  column_to_rownames("pat_id") %>%
  as.matrix()

clinical <- clinical %>%
  filter(pat_id %in% rownames(matrix)) %>%
  arrange(pat_id)

df_metadata <- clinical %>%
  select(cause_of_death, clinical_progr, psa_progr, recurrence, gs_grp, pat_id) %>%
  distinct()

matrix <- matrix[as.character(df_metadata$pat_id), ]

hue_order_list <- list(
  cause_of_death = c("alive", "PCa_death", "non-PCa_death"),
  clinical_progr = c(0, 1),
  psa_progr = c(0, 1),
  recurrence = c("no_recurrence", "recurrence"),
  gs_grp = c("1", "2", "3", "4", "5", "nan"),
  gleason_grp = c("1", "2", "3", "4", "5"),
  inflammation = c("yes", "no", "nan"),
  stromogenic_smc_loss_reactive_stroma_present = c("yes", "no", "nan"),
  glandular_atrophy_pin = c("yes", "no", "nan"),
  ln_status = c(0, 1)
)

for (col in names(hue_order_list)) {
  if (col %in% colnames(df_metadata)) {
    df_metadata[[col]] <- factor(df_metadata[[col]], levels = hue_order_list[[col]])
  }
}

custom_annotation_colors <- list(
  clinical_progr = c("0" = "#457B9D", "1" = "#E63946"),
  cause_of_death = c("alive" = "#06D6A0", "PCa_death" = "#EF476F", "non-PCa_death" = "#457B9D"),
  recurrence = c("no_recurrence" = "#1D3557", "recurrence" = "#F4A261"),
  gleason_grp = c("1" = "yellow", "2" = "#E9C46A", "3" = "#F4A261", "4" = "#E76F51", "5" = "red"),
  inflammation = c("yes" = "#EF476F", "no" = "#457B9D", "nan" = "#999999"),
  stromogenic_smc_loss_reactive_stroma_present = c("yes" = "#E76F51", "no" = "#2A9D8F", "nan" = "#999999"),
  glandular_atrophy_pin = c("yes" = "#F4A261", "no" = "#264653", "nan" = "#999999"),
  ln_status = c("0" = "#2A9D8F", "1" = "#E76F51")
)

annotation_colors <- list()
for (col in colnames(df_metadata)) {
  if (col %in% names(custom_annotation_colors)) {
    annotation_colors[[col]] <- custom_annotation_colors[[col]]
  } else if (col %in% names(hue_order_list)) {
    unique_vals <- hue_order_list[[col]]
    annotation_colors[[col]] <- structure(circlize::rand_color(length(unique_vals)), names = unique_vals)
  }
}

ha <- rowAnnotation(df = df_metadata, col = annotation_colors, show_legend = c(rep(TRUE, 5), FALSE))
col_fun <- colorRamp2(c(0, 1), c("white", "darkgreen"))

jsd_distance <- function(x, y) {
  entropy_x <- entropy(x)
  entropy_y <- entropy(y)
  m <- (x + y) / 2
  entropy(m) - (entropy_x + entropy_y) / 2
}

info_niches <- read.csv(file.path(niches_dir, "annotation", "niche_annotations_v2.csv"))
info_niches <- info_niches %>%
  select(-cluster) %>%
  distinct() %>%
  filter(niche != "unassigned")

matrix <- matrix[, info_niches$niche]

niche_colors <- setNames(unique(info_niches$niche_color), unique(info_niches$niche))
meta_niche_colors <- setNames(unique(info_niches$meta_niche_color), unique(info_niches$meta_niche))

niche_anno <- HeatmapAnnotation(
  df = data.frame(
    meta_niche = factor(info_niches$meta_niche, levels = names(meta_niche_colors)),
    niche = factor(info_niches$niche, levels = names(niche_colors))
  ),
  col = list(meta_niche = meta_niche_colors, niche = niche_colors),
  show_legend = c(TRUE, FALSE)
)

p <- ComplexHeatmap::Heatmap(
  matrix,
  col = col_fun,
  name = "Proportion",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_columns = "spearman",
  clustering_distance_rows = jsd_distance,
  show_row_names = FALSE,
  show_column_names = TRUE,
  right_annotation = ha,
  top_annotation = niche_anno,
  column_split = info_niches$meta_niche,
  row_split = 5,
  heatmap_legend_param = list(title = "Proportion"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8)
)

plot_path <- file.path(save_dir, "figure6a_niche_proportion_heatmap_patient.pdf")
pdf(plot_path, width = 28, height = 20)
ht <- draw(p)
dev.off()
cat("Saved Figure 6a to", plot_path, "\n")

# --- dendrogram-split risk groups + KM (same script in the old repo) ---
row_orders <- row_order(ht)
clustered_rownames <- lapply(row_orders, \(idx) rownames(matrix)[idx])
names(clustered_rownames) <- paste0("cluster_", seq_along(clustered_rownames))
df <- stack(clustered_rownames) %>% rename(pat_id = values, risk_group = ind)

metadata <- merge(df_metadata, df, by = "pat_id", all.x = TRUE)

progression <- clinical %>% select(pat_id, clinical_progr_time) %>% distinct()
death <- clinical %>% select(pat_id, os_status, last_fu) %>% distinct()
death[["overall_survival"]] <- ifelse(death[["os_status"]] == "alive", 0, 1)

metadata <- merge(metadata, progression, by = "pat_id", all.x = TRUE)
metadata <- merge(metadata, death, by = "pat_id", all.x = TRUE)

if (!("risk_group" %in% colnames(metadata))) {
  metadata["risk_group"] <- metadata["cluster"]
}
metadata_filtered <- metadata %>%
  group_by(risk_group) %>%
  filter(n() >= 10) %>%
  ungroup()
metadata_filtered$risk_group <- factor(metadata_filtered$risk_group)
metadata_filtered$clinical_progr <- as.numeric(metadata_filtered$clinical_progr)

fit_prog <- survfit2(Surv(clinical_progr_time, clinical_progr) ~ risk_group, data = metadata_filtered)
p_prog <- fit_prog |>
  ggsurvfit() +
  labs(x = "Time", y = "Progression-free survival probability") +
  add_risktable() +
  add_pvalue()
pdf(file.path(save_dir, "figure6a_progression_km_by_dendrogram_split.pdf"), width = 9, height = 8)
print(p_prog)
dev.off()

fit_os <- survfit2(Surv(last_fu, overall_survival) ~ risk_group, data = metadata_filtered)
p_os <- fit_os |>
  ggsurvfit() +
  labs(x = "Time", y = "Survival probability") +
  add_risktable() +
  add_pvalue()
pdf(file.path(save_dir, "figure6a_survival_km_by_dendrogram_split.pdf"), width = 9, height = 8)
print(p_os)
dev.off()

cat("Saved Figure 6a dendrogram-split KM panels to", save_dir, "\n")
