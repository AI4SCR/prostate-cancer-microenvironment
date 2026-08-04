library(arrow)
result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches"

save_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/kaplan_meier/cell_types/"
dir.create(save_dir, showWarnings = FALSE)

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

df_clusters <- read_parquet("/Users/me3312/Documents/Paper_PCa/5-niches/annotation/clusters_annotated_v2.parquet")
df_clusters[['sample_name']] <- df_clusters[['tma_id']]




######### per niche ###########
df_freqs <- compute_label_frequency(df_clusters, level = "label", pseudocount = 0)

## rename sample_name to tma_id, select niche, proportion and make to wide format
df_wide <- df_freqs %>%
  select(tma_id = sample_name, label, proportion) %>%
  pivot_wider(names_from = label, values_from = proportion, values_fill = 0)

df_props <- df_wide

## plot histogram for selected columns of df_props
cols <- colnames(df_props)[2:ncol(df_props)]
qs <- c(0.25, 0.5, 0.6, 0.75, 0.8, 0.9)

quantiles <- list()
for (col in cols) {
  
  x <- df_props[[col]]
  x_nz <- x[x > 0]
  #x_nz <- x
  
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
  #print(p)
  
  # ggsave(...)
}

clinical.path = file.path('/Users/me3312/Documents/Paper_PCa/0-paper/0-export/clinical.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

tma_ids.valid = intersect(df_props$tma_id, clinical$tma_id)

clinical$os_status <- ifelse(clinical$os_status == "dead", "1", "0")
clinical$disease_progr <- as.character(clinical$disease_progr)
##
clinical_matched = clinical %>%
  filter(tma_id %in% tma_ids.valid) 

##
outcome_cols = c("os_status", "disease_progr")

clinical_tma = clinical_matched %>%
  select(pat_id, all_of(outcome_cols), tma_id) %>%
  distinct()

clinical_patient = clinical_tma %>%
  select(-tma_id) %>%
  group_by(pat_id) %>%
  distinct()
  

outcome_colors = c("0" = "lightblue", "1" = "purple")

df_patient <- clinical_matched %>%
  select(pat_id, tma_id) %>%
  distinct()

freq_patient_max <- df_props %>%
  left_join(df_patient, by = "tma_id") %>%
  select(-tma_id) %>%
  group_by(pat_id) %>%
  summarise(across(where(is.numeric), max, na.rm = TRUE),
            .groups = "drop")

freq_tma <- df_props

###################### on a patient level ############


freq_long <- freq_patient_max %>%
  pivot_longer(cols = -pat_id, names_to = "celltype", values_to = "freq")

clin_long <- clinical_patient %>%
  select(pat_id, all_of(outcome_cols)) %>%
  distinct(pat_id, .keep_all = TRUE) %>%
  pivot_longer(cols = all_of(outcome_cols), names_to = "outcome", values_to = "value") %>%
  mutate(value = as.character(value))
n_pat <- nrow(clinical_patient)


ct_quantiles <- quantiles


for (ct in colnames(df_props)[2:ncol(df_props)]) {
  df_ct <- freq_long %>%
    filter(celltype == ct) %>%
    arrange(freq) %>%
    mutate(pat_pos = row_number())
  
  df_ct <- df_ct %>%
    mutate(freq_plot = ifelse(freq == 0, NA_real_, freq))
  
  clin_ct <- clin_long %>%
    semi_join(df_ct, by = "pat_id") %>%
    left_join(df_ct %>% select(pat_id, pat_pos), by = "pat_id") %>%
    mutate(outcome = factor(outcome, levels = unique(outcome)))
  
  x_breaks <- pretty(seq_len(n_pat), n = 10)
  x_breaks <- x_breaks[x_breaks >= 1 & x_breaks <= n_pat]
  
  # same x-scale for both panels
  x_scale <- scale_x_continuous(
    limits = c(0.5, n_pat + 0.5),
    breaks = x_breaks,
    labels = x_breaks,
    expand = c(0, 0)
  )
  qs <- ct_quantiles[[ct]]
  vlines <- tibble(thresh = qs) %>%
    rowwise() %>%
    mutate(x = if (all(is.na(df_ct$freq))) NA_real_
           else min(df_ct$pat_pos[df_ct$freq >= thresh], na.rm = TRUE)) %>%
    ungroup() %>%
    filter(is.finite(x))
  # ensure vlines are at .5 boundaries if you want them between patients:
  # vlines <- vlines %>% mutate(x = x + 0.5)
  
  p_out <- ggplot(clin_ct, aes(x = pat_pos, y = outcome, fill = value)) +
    geom_tile(height = 0.9) +
    geom_vline(data = vlines, aes(xintercept = x), inherit.aes = FALSE) +   # <--- add here
    x_scale +
    {if (!is.null(outcome_colors)) scale_fill_manual(values = outcome_colors, na.value = "grey90") else NULL} +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x  = element_blank(),     # keep top x labels off
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",           # <--- crucial for equal widths
      plot.margin = margin(2, 6, 0, 6)
    ) +
    labs(title = ct)
  
  p_freq <- ggplot(df_ct, aes(x = pat_pos, y = celltype, fill = freq_plot)) +
    geom_tile(height = 0.9) +
    geom_vline(data = vlines, aes(xintercept = x), inherit.aes = FALSE) +
    x_scale +
    scale_fill_viridis_c(option = "cividis", name = "Frequency", na.value = "grey90") +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x  = element_text(size = 9),
      axis.ticks.x = element_line(),
      panel.grid = element_blank(),
      legend.position = "right",           # <--- crucial for equal widths
      plot.margin = margin(0, 6, 2, 6)
    )
  
  print(p_out / p_freq + plot_layout(heights = c(1, 0.6)))
  p_save <- p_out / p_freq + plot_layout(heights = c(1, 0.6))
  plot_dir = file.path(save_dir, "diagrams" , "patient_level")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = file.path(plot_dir, paste0("kaplan_meier_", ct, "_pat_level.png")),
         plot = p_save,
         width = 14,
         height = 12,
         dpi = 300)

  
}


###################### on a tma level ############


freq_long <- freq_tma %>%
  pivot_longer(cols = -tma_id, names_to = "celltype", values_to = "freq")

clin_long <- clinical_tma %>%
  select(tma_id, all_of(outcome_cols)) %>%
  distinct(tma_id, .keep_all = TRUE) %>%
  pivot_longer(cols = all_of(outcome_cols), names_to = "outcome", values_to = "value") %>%
  mutate(value = as.character(value))

n_pat <- nrow(clinical_tma)

ct_quantiles <- quantiles

for (ct in colnames(df_props)[2:ncol(df_props)]) {
  df_ct <- freq_long %>%
    filter(celltype == ct) %>%
    arrange(freq) %>%
    mutate(pat_pos = row_number())
  df_ct <- df_ct %>%
    mutate(freq_plot = ifelse(freq == 0, NA_real_, freq))
  
  clin_ct <- clin_long %>%
    semi_join(df_ct, by = "tma_id") %>%
    left_join(df_ct %>% select(tma_id, pat_pos), by = "tma_id") %>%
    mutate(outcome = factor(outcome, levels = unique(outcome)))
  
  x_breaks <- pretty(seq_len(n_pat), n = 10)
  x_breaks <- x_breaks[x_breaks >= 1 & x_breaks <= n_pat]
  
  # same x-scale for both panels
  x_scale <- scale_x_continuous(
    limits = c(0.5, n_pat + 0.5),
    breaks = x_breaks,
    labels = x_breaks,
    expand = c(0, 0)
  )
  qs <- ct_quantiles[[ct]]
  vlines <- tibble(thresh = qs) %>%
    rowwise() %>%
    mutate(x = if (all(is.na(df_ct$freq))) NA_real_
           else min(df_ct$pat_pos[df_ct$freq >= thresh], na.rm = TRUE)) %>%
    ungroup() %>%
    filter(is.finite(x))
  # ensure vlines are at .5 boundaries if you want them between patients:
  # vlines <- vlines %>% mutate(x = x + 0.5)
  
  p_out <- ggplot(clin_ct, aes(x = pat_pos, y = outcome, fill = value)) +
    geom_tile(height = 0.9) +
    geom_vline(data = vlines, aes(xintercept = x), inherit.aes = FALSE) +   # <--- add here
    x_scale +
    {if (!is.null(outcome_colors)) scale_fill_manual(values = outcome_colors, na.value = "grey90") else NULL} +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x  = element_blank(),     # keep top x labels off
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",           # <--- crucial for equal widths
      plot.margin = margin(2, 6, 0, 6)
    ) +
    labs(title = ct)
  
  p_freq <- ggplot(df_ct, aes(x = pat_pos, y = celltype, fill = freq_plot)) +
    geom_tile(height = 0.9) +
    geom_vline(data = vlines, aes(xintercept = x), inherit.aes = FALSE) +
    x_scale +
    scale_fill_viridis_c(option = "cividis", name = "Frequency", na.value = "grey90") +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x  = element_text(size = 9),
      axis.ticks.x = element_line(),
      panel.grid = element_blank(),
      legend.position = "right",           # <--- crucial for equal widths
      plot.margin = margin(0, 6, 2, 6)
    )
  
  print(p_out / p_freq + plot_layout(heights = c(1, 0.6)))
  # save plot
  p_save <- p_out / p_freq + plot_layout(heights = c(1, 0.6))
  plot_dir = file.path(save_dir, "diagrams")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = file.path(plot_dir, paste0("kaplan_meier_", ct, "_tma_level.png")),
         plot = p_save,
         width = 14,
         height = 12,
         dpi = 300)
  
  
  
}





