# Reviewer response (Issue 2, cell-type level): does niche 9
# (luminal_CAF1(CD105High)) remain associated with outcome once the
# proportion of its defining stromal cell type, stromal-CAF1(CD105+), is in
# the model? Both terms are CLR-transformed proportions on the same
# tumor-cores-only, max-pooled-per-patient pipeline as
# figure4de_cox_hazard_ratio.R -- one built from `cell_annotation.parquet`'s
# `niche` column, the other from `metadata.parquet`'s `label` column.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(compositions)
library(ggplot2)
library(tibble)
library(readr)

base_dir <- Sys.getenv("BASE_DIR")
export_dir <- Sys.getenv("EXPORT_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(export_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

save_dir <- file.path(dirname(output_figures_dir), "revision", "figure7_niche9")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

niche_col <- "luminal_CAF1(CD105High)" # niche 9
cell_type_col <- "stromal-CAF1(CD105+)" # CAF1 CD105-high

clinical_full <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical_full$pat_id |> n_distinct()

clinical.names <- c("sample_id", "pat_id", "tma_id", "last_fu", "os_status", "disease_progr", "disease_progr_time")
clinical <- clinical_full |> select(any_of(clinical.names))

sample_tma <- clinical_full |>
  select(sample_id, tma_id, pat_id, is_tumor) |>
  distinct()
tma_pat <- sample_tma |> select(tma_id, pat_id) |> distinct()

clr_composition <- function(cells, level_col) {
  # tumor-cores-only, pseudocount-1 freq per tma_id, max-pooled to patient,
  # CLR across the full composition -- mirrors figure4de_cox_hazard_ratio.R
  freq <- cells |>
    count(tma_id, .data[[level_col]], name = "n") |>
    complete(tma_id, .data[[level_col]], fill = list(n = 0)) |>
    mutate(n = n + 1) |>
    group_by(tma_id) |>
    mutate(prop = n / sum(n)) |>
    ungroup()

  composition <- freq |>
    select(tma_id, all_of(level_col), prop) |>
    left_join(tma_pat, by = "tma_id") |>
    group_by(pat_id, .data[[level_col]]) |>
    summarise(score = max(prop), .groups = "drop") |>
    pivot_wider(names_from = all_of(level_col), values_from = score)

  names <- setdiff(colnames(composition), "pat_id")
  composition[, names] <- clr(as.matrix(composition[, names]))
  composition
}

# %% niche composition (from cell_annotation.parquet)
cells_niche <- read_parquet(file.path(export_dir, "cell_annotation.parquet")) |>
  select(sample_id, object_id, niche) |>
  inner_join(sample_tma, by = "sample_id") |>
  filter(!is.na(is_tumor), is_tumor == "yes")
niche_composition <- clr_composition(cells_niche, "niche")
stopifnot("niche 9 column not found in niche composition" = niche_col %in% colnames(niche_composition))
niche_data <- niche_composition |> select(pat_id, niche9_clr = all_of(niche_col))

# %% cell-type composition (from metadata.parquet)
metadata <- read_parquet(file.path(export_dir, "metadata.parquet"))
cells_label <- metadata |>
  select(sample_id, object_id, label) |>
  inner_join(sample_tma, by = "sample_id") |>
  filter(!is.na(is_tumor), is_tumor == "yes")
label_composition <- clr_composition(cells_label, "label")
stopifnot("CAF1(CD105+) column not found in cell-type composition" = cell_type_col %in% colnames(label_composition))
cell_type_data <- label_composition |> select(pat_id, caf1cd105_clr = all_of(cell_type_col))

cox_data_base <- niche_data |> inner_join(cell_type_data, by = "pat_id")

fit_cox <- function(df, formula_rhs) {
  fml <- as.formula(paste0("Surv(time, event) ~ ", formula_rhs))
  fit <- coxph(fml, data = df)
  s <- summary(fit)
  terms <- rownames(s$coefficients)
  tibble(
    term = terms,
    hr = unname(s$coefficients[, "exp(coef)"]),
    hr_lower = s$conf.int[, "lower .95"],
    hr_upper = s$conf.int[, "upper .95"],
    p_value = unname(s$coefficients[, "Pr(>|z|)"]),
    n = s$n,
    n_event = s$nevent
  )
}

for (spec in list(
  list(event.name = "os_status", label = "os"),
  list(event.name = "disease_progr", label = "progr")
)) {
  event.name <- spec$event.name
  label <- spec$label

  if (event.name == "os_status") {
    event.value <- "dead"
    time.name <- "last_fu"
  } else {
    event.value <- 1
    time.name <- "disease_progr_time"
  }

  clinical.surv <- clinical |>
    select(pat_id, all_of(c(event.name, time.name))) |>
    distinct()
  stopifnot(nrow(clinical.surv) == num.patients)
  clinical.surv["event"] <- clinical.surv[[event.name]] == event.value
  clinical.surv["time"] <- clinical.surv[[time.name]]

  cox_data <- cox_data_base |> inner_join(clinical.surv, by = "pat_id")
  stopifnot(!(is.na(cox_data$event) |> any()))
  stopifnot(!(is.na(cox_data$time) |> any()))

  univariate_niche <- fit_cox(cox_data, "niche9_clr") |> mutate(model = "niche9_only")
  univariate_cell <- fit_cox(cox_data, "caf1cd105_clr") |> mutate(model = "caf1cd105_only")
  multivariate <- fit_cox(cox_data, "niche9_clr + caf1cd105_clr") |> mutate(model = "niche9_plus_caf1cd105")

  results <- bind_rows(univariate_niche, univariate_cell, multivariate)
  results$outcome <- event.name

  save.path <- file.path(save_dir, sprintf("figure7_cox_niche9_caf1cd105_%s.csv", label))
  write_csv(results, save.path)

  pdat <- multivariate |>
    mutate(term = factor(term, levels = term))

  p <- ggplot(pdat, aes(y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), height = 0.2, linewidth = 0.6, color = "black") +
    geom_point(aes(x = hr), size = 2.5, shape = 21, fill = "red") +
    labs(
      x = "Hazard ratio",
      y = NULL,
      title = paste("Niche 9 + CAF1(CD105+) -- multivariate Cox --", event.name),
      subtitle = "Points = HR; bars = 95% CI"
    ) +
    theme_minimal(base_size = 11)
  ggsave(
    filename = file.path(save_dir, sprintf("figure7_cox_niche9_caf1cd105_%s.png", label)),
    plot = p, width = 6, height = 4
  )
}

cat("Saved niche 9 + CAF1(CD105+) Cox results to", save_dir, "\n")
