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

figures_dir <- file.path(output_figures_dir, "figure6")
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
df_cells <- read_parquet(file.path(export_dir, "metadata.parquet"))

sample_col <- "tma_id"
df_cells[["sample_name"]] <- df_cells[[sample_col]]

# compute label frequencies for each sample and label
label_col <- "niche"
df_freqs <- compute_label_frequency(df_cells, level = label_col, pseudocount = 0)
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
plot_name <- paste0("kaplan_meier_inflammation_os_", "risk_group_max.pdf")
ggsave(filename = file.path(figures_dir, plot_name), plot = p1$plot, width = 8, height = 6, dpi = 300)
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
plot_name <- paste0("kaplan_meier_inflammation_progression_", "risk_group_max.pdf")
ggsave(filename = file.path(figures_dir, plot_name), plot = p2$plot, width = 8, height = 6, dpi = 300)
