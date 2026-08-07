# Reviewer response (Issue 2): does niche 9 (luminal_CAF1(CD105High)) remain
# associated with outcome once patient-level Gleason grade (gs_grp) is in the
# model, vs a Gleason-grade-only model? Extends figure7b_niche9_km.R (KM/binary
# risk-group) and reuses figure4de_cox_hazard_ratio.R's CLR-abundance pipeline,
# applied to niches instead of cell types, so niche 9 enters as a continuous
# covariate rather than a median-split group.

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
data_dir <- Sys.getenv("DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(data_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot(
  "refusing to treat BASE_DIR as writable" = !startsWith(normalizePath(data_dir, mustWork = FALSE), normalizePath(base_dir, mustWork = FALSE))
)

save_dir <- file.path(dirname(output_figures_dir), "revision", "figure7_niche9")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

niche_col <- "luminal_CAF1(CD105High)" # niche 9

clinical_full <- read_parquet(file.path(data_dir, "clinical.parquet"))
num.patients <- clinical_full$pat_id |> n_distinct()

clinical.names <- c("sample_id", "pat_id", "tma_id", "last_fu", "os_status", "disease_progr", "disease_progr_time", "gs_grp")
clinical <- clinical_full |> select(any_of(clinical.names))

# %% per-tma_id niche composition, restricted to tumor cores, max-pooled per
# patient -- same pipeline as figure4de_cox_hazard_ratio.R applied to cell
# types, here applied to `cell_annotation.parquet`'s `niche` column instead of
# `metadata.parquet`'s `label` column
cells_niche <- read_parquet(file.path(data_dir, "cells", "cell_annotation.parquet")) |>
  select(sample_id, object_id, niche)
sample_tma <- clinical_full |>
  select(sample_id, tma_id, pat_id, is_tumor) |>
  distinct()

cells <- cells_niche |>
  inner_join(sample_tma, by = "sample_id") |>
  filter(!is.na(is_tumor), is_tumor == "yes")

freq <- cells |>
  count(tma_id, niche, name = "n") |>
  complete(tma_id, niche, fill = list(n = 0)) |>
  mutate(n = n + 1) |>
  group_by(tma_id) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

tma_pat <- sample_tma |> select(tma_id, pat_id) |> distinct()

composition <- freq |>
  select(tma_id, niche, prop) |>
  left_join(tma_pat, by = "tma_id") |>
  group_by(pat_id, niche) |>
  summarise(score = max(prop), .groups = "drop") |>
  pivot_wider(names_from = niche, values_from = score)

niche_names <- setdiff(colnames(composition), "pat_id")
stopifnot("niche 9 column not found in niche composition" = niche_col %in% niche_names)

data <- composition
data[, niche_names] <- clr(as.matrix(data[, niche_names]))
data <- data |> select(pat_id, niche9_clr = all_of(niche_col))

# %% patient-level Gleason grade group (gs_grp), ordinal 1-5
gleason <- clinical_full |>
  select(pat_id, gs_grp) |>
  distinct()
stopifnot("gs_grp is not one value per patient" = nrow(gleason) == num.patients)
gleason$gs_grp <- suppressWarnings(as.numeric(gleason$gs_grp))
n_gleason_na <- sum(is.na(gleason$gs_grp))
cat("Patients with missing gs_grp (dropped from these models):", n_gleason_na, "\n")

cox_data_base <- data |>
  inner_join(gleason, by = "pat_id") |>
  filter(!is.na(gs_grp))

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

  univariate <- fit_cox(cox_data, "gs_grp") |> mutate(model = "gleason_only")
  multivariate <- fit_cox(cox_data, "niche9_clr + gs_grp") |> mutate(model = "niche9_plus_gleason")

  results <- bind_rows(univariate, multivariate)
  results$outcome <- event.name
  results$n_patients_gs_grp_dropped <- n_gleason_na

  save.path <- file.path(save_dir, sprintf("figure7_cox_niche9_gleason_%s.csv", label))
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
      title = paste("Niche 9 + Gleason grade -- multivariate Cox --", event.name),
      subtitle = "Points = HR; bars = 95% CI"
    ) +
    theme_minimal(base_size = 11)
  ggsave(
    filename = file.path(save_dir, sprintf("figure7_cox_niche9_gleason_%s.png", label)),
    plot = p, width = 6, height = 4
  )
}

cat("Saved niche 9 + Gleason Cox results to", save_dir, "\n")
