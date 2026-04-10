library(arrow)
result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches"

base_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies"
df_props = read_parquet(file.path(base_dir, "stacked_barplots/props_niche_tma_id.parquet"))
# metadata = read_parquet(file.path(base_dir, "metadata_aligned_niche_frequencies.parquet"))
# df_props = read_parquet(file.path(result_dir, "frequencies/stacked_barplots/props_tma.parquet"))
#clusters = read_parquet(file.path(result_dir, "frequencies/stacked_barplots/clustered_metadata.parquet"))
# ln_status = read_parquet(file.path(result_dir, "frequencies/ln_status.parquet"))
# clinical = read_parquet(file.path(result_dir, "frequencies/clinical.parquet"))
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(entropy)
library(survival)
library(survminer)

clinical.path = file.path('/Users/me3312/Documents/Paper_PCa/0-paper/0-export/clinical.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

tma_ids.valid = intersect(df_props$tma_id, clinical$tma_id)



compute_label_frequency <- function(data, level, pseudocount = 1) {
  level_sym <- rlang::sym(level)
  
  data_summary <- data %>%
    dplyr::group_by(sample_name, !!level_sym) %>%
    dplyr::summarise(count = dplyr::n()+pseudocount, .groups = "drop") %>%
    tidyr::complete(
      sample_name,
      !!level_sym,
      fill = list(count = 0)
    )
  
  df_freqs <- data_summary %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(
      proportion = (count) / (sum(count))
    ) %>%
    dplyr::ungroup()
  
  return(df_freqs)
}

df_clusters <- read_parquet("/Users/me3312/Documents/Paper_PCa/5-niches/annotation/clusters_annotated_v2.parquet")
df_clusters[['sample_name']] <- df_clusters[['tma_id']]






######### per niche ###########
df_freqs <- compute_label_frequency(df_clusters, level = "niche", pseudocount = 0)

## rename sample_name to tma_id, select niche, proportion and make to wide format
df_wide <- df_freqs %>%
  select(tma_id = sample_name, niche, proportion) %>%
  pivot_wider(names_from = niche, values_from = proportion, values_fill = 0)

df_props <- df_wide

## plot histogram for selected columns of df_props
cols_inflamed <- c("TLS", "Macrophages_Tcells_CAF1(CD105-)", "immune_bloodvessels_CAF1(CD105-)")
qs <- c(0.25, 0.5, 0.6, 0.75, 0.8)

quantiles <- list()

for (col in cols_inflamed) {
  
  x <- df_props[[col]]
  x_nz <- x[x > 0]
  
  q_vals <- quantile(x_nz, probs = qs, na.rm = TRUE)
  print(paste("Quantiles for", col, ":"))
  print(q_vals)
  quantiles[[col]] <- q_vals
  
  p <- ggplot(df_props, aes(x = .data[[col]])) +
    geom_histogram(
      binwidth = 0.01,
      fill = "blue",
      color = "black",
      alpha = 0.7
    ) +
    geom_vline(
      xintercept = q_vals,
      linetype = "dashed",
      linewidth = 1,
      color = "red"
    ) +
    labs(
      title = paste("Histogram of", col),
      x = col,
      y = "Frequency"
    ) +
    theme_minimal()
  print(p)
  
  # ggsave(...)
}

threshold <- "75%"
thr <- sapply(cols_inflamed, function(col) quantiles[[col]][[threshold]])
names(thr) <- cols_inflamed

df_inflammation <- df_props %>%
  select(tma_id, all_of(cols_inflamed)) %>%
  mutate(
    n_inflamed = rowSums(
      sweep(across(all_of(cols_inflamed)), 2, thr, `>`),
      na.rm = TRUE
    ),
    inflammation = ifelse(n_inflamed > 1, 1, 0)
  ) %>%
  select(tma_id, n_inflamed, inflammation)

# join with clinical to get patient id


clinical.path = file.path('/Users/me3312/Documents/Paper_PCa/0-paper/0-export/clinical.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

tma_ids.valid = intersect(df_props$tma_id, clinical$tma_id)

clinical$os_status <- ifelse(clinical$os_status == "dead", 1, 0)

progression <- clinical %>%
  select(
    pat_id,
    disease_progr,
    disease_progr_time,
  ) %>%
  distinct()

death <- clinical %>%
  select(
    pat_id,
    os_status,
    last_fu
  ) %>%
  distinct()

df_patient <- clinical %>%
  select(
    pat_id,
    tma_id
  ) %>%
  distinct() %>%
  inner_join(df_inflammation, by = "tma_id") %>%
  select(pat_id, risk_group = n_inflamed)

df_patient <- df_patient %>%
  group_by(pat_id) %>%
  summarise(
    risk_group = max(risk_group),
    .groups = "drop"
  )

df_analysis <- df_patient %>%
  inner_join(progression, by = "pat_id") %>%
  inner_join(death, by = "pat_id")

save_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/kaplan_meier/inflammation/"
dir.create(save_dir, showWarnings = FALSE)


fit <- survfit(Surv(last_fu, os_status) ~ risk_group, data = df_analysis)
p1 <- ggsurvplot(
  fit,
  data = df_analysis,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = "Set2",
  xlab = "Time",
  ylab = "Survival probability using 75% quantile",
  legend.title = "Group",
  risk.table.height = 0.25
)
plot_name <- paste0("kaplan_meier_inflammation_os_", "risk_group_bin.pdf")
#ggsave(filename = file.path(save_dir, plot_name), plot = p1$plot, width = 8, height = 6, dpi = 300)
print(p1)
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
  risk.table.height = 0.25
)
plot_name <- paste0("kaplan_meier_inflammation_progression_", "risk_group_bin.pdf")
#ggsave(filename = file.path(save_dir, plot_name), plot = p2$plot, width = 8, height = 6, dpi = 300)
print(p2)

df_histo <- clinical %>%
  select(
    tma_id,
    inflammation,
    sample_name
  ) %>%
  rename(eva_annotation_infl = inflammation) %>%
  distinct() %>%
  inner_join(df_inflammation, by = "tma_id")

## save as csv
write.csv(df_histo, file.path(save_dir, "histo_inflammation_annotation.csv"), row.names = FALSE)
