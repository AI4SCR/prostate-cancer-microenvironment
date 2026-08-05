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
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))

results_dir <- file.path(output_figures_dir, "figureS6", "niches")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
figures_dir <- file.path(output_figures_dir, "figureS6", "niches")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

compute_label_frequency <- function(data, level, pseudocount = 1) {
  level_sym <- rlang::sym(level)

  data_summary <- data %>%
    dplyr::group_by(sample_name, !!level_sym) %>%
    dplyr::summarise(count = dplyr::n() + pseudocount, .groups = "drop") %>%
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

### read clinical metadata
clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()
print(paste("Number of patients:", num.patients))

### read cell annotation data
df_cells <- read_parquet(file.path(export_dir, "cell_annotation.parquet"))

sample_col <- "tma_id"
df_cells[["sample_name"]] <- df_cells[[sample_col]]

# compute label frequencies for each sample and label
label_col <- "niche"
df_freqs <- compute_label_frequency(df_cells, level = label_col, pseudocount = 0)
df_wide <- df_freqs %>%
  select(tma_id = sample_name, niche, proportion) %>%
  pivot_wider(names_from = niche, values_from = proportion, values_fill = 0)

df_props <- df_wide

## plot histogram for celltype columns of df_props
cols <- c("luminal_infiltrated", "tumor_CAF1_lymphocytes") # niches 2 and 8
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
  # print(p)
}

## set median as threshold for binary classification
threshold <- "50%"
thr <- sapply(cols, function(col) quantiles[[col]][[threshold]])
names(thr) <- cols

# create binary labels based on threshold at sample level
df_binary <- df_props %>%
  mutate(across(all_of(cols), ~ ifelse(.x >= thr[cur_column()], 1, 0)))

# join with clinical to get patient id and survival data
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

### run survival analysis for each label and save results
for (col in cols) {
  df_label <- df_binary
  df_label[["target"]] <- df_label[[col]]
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
  form <- as.formula(fit$call$formula)
  sd <- survdiff(form, data = df_analysis)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  cox <- coxph(formula = form, data = df_analysis)
  s_cox <- summary(cox)$coefficients[, c("exp(coef)", "Pr(>|z|)")]
  results_os[[col]] <- list(pval = pval, cox_coef = s_cox[1], cox_pval = s_cox[2])
  p_os <- ggsurvplot(
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
  # print(p_os)

  fit_prog <- survfit(Surv(disease_progr_time, disease_progr) ~ risk_group, data = df_analysis)
  form <- as.formula(fit_prog$call$formula)
  sd <- survdiff(form, data = df_analysis)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  cox <- coxph(formula = form, data = df_analysis)
  s_cox <- summary(cox)$coefficients[, c("exp(coef)", "Pr(>|z|)")]
  results_prog[[col]] <- list(pval = pval, cox_coef = s_cox[1], cox_pval = s_cox[2])
  p_prog <- ggsurvplot(
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
  # print(p_prog)
  if (col %in% c("luminal_infiltrated", "luminal_CAF1(CD105High)", "tumor_CAF1(CD105High)")) {
    p_prog_combined <- p_prog$plot / p_prog$table
    plot_name <- paste0(figures_dir, "km_survival_", "_disease_progr_", col, "with_table.pdf")
    # ggsave(plot_name, p_prog_combined, width = 8, height = 6, dpi = 300)
  }
}

## aggregate results into dataframes and adjust p-values for multiple testing (BH method)
df_os <- do.call(rbind, lapply(names(results_os), function(col) {
  data.frame(
    niche = col,
    pval = results_os[[col]]$pval,
    cox_coef = results_os[[col]]$cox_coef,
    cox_pval = results_os[[col]]$cox_pval
  )
}))
df_os$qval <- round(p.adjust(df_os$pval, method = "BH"), 4)
df_os$cox_qval <- round(p.adjust(df_os$cox_pval, method = "BH"), 4)

df_prog <- do.call(rbind, lapply(names(results_prog), function(col) {
  data.frame(
    niche = col,
    pval = results_prog[[col]]$pval,
    cox_coef = results_prog[[col]]$cox_coef,
    cox_pval = results_prog[[col]]$cox_pval
  )
}))
df_prog$qval <- round(p.adjust(df_prog$pval, method = "BH"), 4)
df_prog$cox_qval <- round(p.adjust(df_prog$cox_pval, method = "BH"), 4)

write.csv(df_os, file.path(results_dir, "overall_survival_analysis_results.csv"), row.names = FALSE)
write.csv(df_prog, file.path(results_dir, "progression_free_survival_analysis_results.csv"), row.names = FALSE)
