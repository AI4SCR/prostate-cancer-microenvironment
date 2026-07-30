library(arrow)
result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches"

save_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/kaplan_meier/niches/"
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
df_freqs <- compute_label_frequency(df_clusters, level = "niche", pseudocount = 0)

## rename sample_name to tma_id, select niche, proportion and make to wide format
df_wide <- df_freqs %>%
  select(tma_id = sample_name, niche, proportion) %>%
  pivot_wider(names_from = niche, values_from = proportion, values_fill = 0)

df_props <- df_wide

## plot histogram for selected columns of df_props
cols <- colnames(df_props)[2:ncol(df_props)]
qs <- c(0.25, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)

quantiles <- list()

for (col in cols) {
  
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

threshold <- "50%"
thr <- sapply(cols, function(col) quantiles[[col]][[threshold]])
names(thr) <- cols


save_dir = paste(save_dir, "threshold_", substr(threshold, 1, 2), "/", sep = "")
dir.create(save_dir, showWarnings = FALSE)

progression_dir = paste0(save_dir, "progression_free_survival/")
dir.create(progression_dir, showWarnings = FALSE)
survival_dir = paste0(save_dir, "overall_survival/")
dir.create(survival_dir, showWarnings = FALSE)

#
df_binary <- df_props %>%
  mutate(across(-tma_id, ~ ifelse(.x >= thr[cur_column()], 1, 0)))

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


results_os <- list()
results_prog <- list()

for (col in cols){
  df_niche <- df_binary 
  df_niche[['target']] <- df_niche[[col]]
  df_patient <- clinical %>%
    select(
      pat_id,
      tma_id
    ) %>%
    distinct() %>%
    inner_join(df_niche, by = "tma_id") %>%
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
  form <- as.formula(fit$call$formula)
  sd <- survdiff(form, data = df_analysis)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  cox <- coxph(formula = form, data = df_analysis)
  s_cox <- summary(cox)$coefficients[, c("exp(coef)", "Pr(>|z|)")]
  results_os[[col]] <- list(pval = pval, cox_coef = s_cox[1], cox_pval = s_cox[2])
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
  plot_name = paste0(survival_dir, "km_survival_", "os_status_", col, ".pdf")
  #ggsave(plot_name, p1$plot, width = 8, height = 6, dpi = 300)
  
  fit_prog <- survfit(Surv(disease_progr_time, disease_progr) ~ risk_group, data = df_analysis)
  form <- as.formula(fit_prog$call$formula)
  sd <- survdiff(form, data = df_analysis)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  cox <- coxph(formula = form, data = df_analysis)
  s_cox <- summary(cox)$coefficients[, c("exp(coef)", "Pr(>|z|)")]
  results_prog[[col]] <- list(pval = pval, cox_coef = s_cox[1], cox_pval = s_cox[2])
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
  plot_name = paste0(progression_dir, "km_progression_", "disease_progr_", col, ".pdf")
  #ggsave(plot_name, p2$plot, width = 8, height = 6, dpi = 300)
  print(p2)
}

df_os = do.call(rbind, lapply(names(results_os), function(col) {
  data.frame(
    niche = col,
    pval = results_os[[col]]$pval,
    cox_coef = results_os[[col]]$cox_coef,
    cox_pval = results_os[[col]]$cox_pval
  )
}))
df_os$qval <- round(p.adjust(df_os$pval, method = "BH"), 4)
df_os$cox_qval <- round(p.adjust(df_os$cox_pval, method = "BH"),4)


df_prog = do.call(rbind, lapply(names(results_prog), function(col) {
  data.frame(
    niche = col,
    pval = results_prog[[col]]$pval,
    cox_coef = results_prog[[col]]$cox_coef,
    cox_pval = results_prog[[col]]$cox_pval
  )
}))
df_prog$qval <- round(p.adjust(df_prog$pval, method = "BH"), 4)
df_prog$cox_qval <- round(p.adjust(df_prog$cox_pval, method = "BH"),4)
