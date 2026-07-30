library(arrow)
# result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches"

save_dir = "/Users/adrianomartinelli/Downloads/me3312"
dir.create(save_dir, showWarnings = FALSE)

base_dir = "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export"
# df_props = read_parquet(file.path(base_dir, "stacked_barplots/props_niche_tma_id.parquet"))
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
# library(entropy)
library(survival)
library(survminer)

clinical.path = file.path(base_dir, 'clinical.parquet')
metadata.path = file.path(base_dir, 'metadata.parquet')

clinical = read_parquet(clinical.path)
metadata = read_parquet(metadata.path)
num.patients = clinical$pat_id |> n_distinct()

metadata <- metadata %>%
  left_join(
    clinical %>% select(sample_id, tma_id, is_tumor),
    by = "sample_id"
  )
stopifnot(metadata$sample_id |> is.na() |> sum() == 0)
stopifnot(metadata$tma_id |> is.na() |> sum() == 0)

tma_ids.valid = intersect(metadata$tma_id, clinical$tma_id)
length(tma_ids.valid)

tumors_only = TRUE
if(tumors_only){
  metadata = metadata |> filter(is_tumor == 'yes')
}
stopifnot(all(!is.na(metadata$is_tumor)))
table(metadata$is_tumor)

compute_label_frequency <- function(data, level, pseudocount = 1) {
  level_sym <- rlang::sym(level)
  
  data_summary <- data %>%
    dplyr::group_by(sample_name, !!level_sym) %>%
    dplyr::summarise(count = dplyr::n()+pseudocount, .groups = "drop") %>%
    tidyr::complete(
      sample_name,
      !!level_sym,
      fill = list(count = pseudocount)
    )
  
  df_freqs <- data_summary %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(
      proportion = (count) / (sum(count))
    ) %>%
    dplyr::ungroup()
  
  return(df_freqs)
}

# df_clusters <- read_parquet("/Users/me3312/Documents/Paper_PCa/5-niches/annotation/clusters_annotated_v2.parquet")
metadata['sample_name'] <- metadata$tma_id


######### per niche ###########
# NOTE: this is with pseudocount == 0 not 1
df_freqs <- compute_label_frequency(metadata, level = "label", pseudocount = 0)

## rename sample_name to tma_id, select niche, proportion and make to wide format
df_wide <- df_freqs %>%
  select(tma_id = sample_name, label, proportion) %>%
  pivot_wider(names_from = label, values_from = proportion, values_fill = 0)

df_props <- df_wide
stopifnot(nrow(df_props) == 459)

# df_props = df_props[order(df_props$tma_id), order(colnames(df_props))]
# head(df_props)
# write_csv(df_props, file.path(save_dir, "proportions-melissa.csv"))

## plot histogram for selected columns of df_props
cols <- colnames(df_props)[2:ncol(df_props)]
qs <- c(0.25, 0.5, 0.6, 0.75, 0.8)

quantiles <- list()

for (col in cols) {
  
  x <- df_props[[col]]
  x_nz <- x[x > 0]
  x_nz <- x
  
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
  # print(p)
  
  # ggsave(...)
}

threshold <- "75%"
thr <- sapply(cols, function(col) quantiles[[col]][[threshold]])
names(thr) <- cols

#
df_binary <- df_props %>%
  mutate(across(-tma_id, ~ ifelse(.x >= thr[cur_column()], 1, 0)))

# join with clinical to get patient id

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

col = "stromal-CAF1(CD105+)"
cols = c(col)
for (col in cols){
  df_label <- df_binary 
  df_label[['target']] <- df_label[[col]]
  df_patient <- clinical %>%
    select(
      pat_id,
      tma_id
    ) %>%
    distinct() %>%
    inner_join(df_label, by = "tma_id") %>%
    select(pat_id, risk_group = target)
  
  df_patient <- df_patient %>%
    group_by(pat_id) %>%
    summarise(
      risk_group = max(risk_group),
      .groups = "drop"
    )
  
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
  print(p1)
  plot_name = file.path(save_dir, paste0("km_survival_", "_os_status_", col, ".pdf"))
  ggsave(plot_name, p1$plot, width = 8, height = 6, dpi = 300)
  
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
  print(p2)
  plot_name = file.path(save_dir, paste0("km_survival_", "_disease_progr_", col, ".pdf"))
  ggsave(plot_name, p2$plot, width = 8, height = 6, dpi = 300)
}
