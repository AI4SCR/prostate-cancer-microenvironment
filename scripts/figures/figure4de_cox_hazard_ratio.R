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

save_dir <- file.path(output_figures_dir, "figure4")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

score_type <- "proportion_tma"
aggregation <- "max"
tumors_only <- TRUE

clinical_full <- read_parquet(file.path(data_dir, "clinical.parquet"))
num.patients <- clinical_full$pat_id |> n_distinct()

clinical.names <- c("sample_id", "pat_id", "tma_id", "last_fu", "os_status", "disease_progr", "disease_progr_time")
clinical <- clinical_full |> select(any_of(clinical.names))

# %% per-tma_id cell-type composition, restricted to tumor cores (tumors_only=TRUE),
# then max-pooled per patient (aggregation='max') -- see docstring
metadata <- read_parquet(file.path(data_dir, "cells", "metadata.parquet"))
sample_tma <- clinical_full |>
  select(sample_id, tma_id, pat_id, is_tumor) |>
  distinct()

cells <- metadata |>
  select(sample_id, object_id, label) |>
  inner_join(sample_tma, by = "sample_id") |>
  filter(!is.na(is_tumor), is_tumor == "yes")

# pseudocount=1 frequency table per tma_id (matches figure4_patient_clustering.py's
# compute_label_frequency, there grouped by pat_id, here by tma_id)
freq <- cells |>
  count(tma_id, label, name = "n") |>
  complete(tma_id, label, fill = list(n = 0)) |>
  mutate(n = n + 1) |>
  group_by(tma_id) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

tma_pat <- sample_tma |> select(tma_id, pat_id) |> distinct()

# max-pool per patient across their tma_ids (aggregation == 'max')
composition <- freq |>
  select(tma_id, label, prop) |>
  left_join(tma_pat, by = "tma_id") |>
  group_by(pat_id, label) |>
  summarise(score = max(prop), .groups = "drop") |>
  pivot_wider(names_from = label, values_from = score)

score_names <- setdiff(colnames(composition), "pat_id")

data <- composition
data[, score_names] <- clr(as.matrix(data[, score_names]))

cox.fit <- function(data, score_name) {
  cols <- c("time", "event", score_name)
  df <- data |> select(all_of(cols))

  fml <- as.formula(paste0("Surv(time, event) ~ `", score_name, "`"))
  fit <- coxph(fml, data = df)
  s <- summary(fit)

  hr <- unname(s$coefficients[1, "exp(coef)"])
  hr.lower <- s$conf.int[1, "lower .95"]
  hr.upper <- s$conf.int[1, "upper .95"]
  p_value <- unname(s$coefficients[1, "Pr(>|z|)"])

  tibble(
    score_name = score_name,
    hr = hr,
    hr_lower = hr.lower,
    hr_upper = hr.upper,
    p_value = p_value,
  )
}

# event.name toggle -- looped over both values instead of a manual re-run (see docstring)
for (spec in list(
  list(event.name = "os_status", label = "e"),
  list(event.name = "disease_progr", label = "d")
)) {
  event.name <- spec$event.name
  label <- spec$label

  if (event.name == "os_status") {
    event.value <- "dead"
    time.name <- "last_fu"
  } else if (event.name == "disease_progr") {
    event.value <- 1
    time.name <- "disease_progr_time"
  }

  cols.surv <- c("pat_id", event.name, time.name)
  clinical.surv <- clinical |> select(all_of(cols.surv)) |> distinct()
  stopifnot(nrow(clinical.surv) == num.patients)

  clinical.surv["event"] <- clinical.surv[[event.name]] == event.value
  clinical.surv["time"] <- clinical.surv[[time.name]]

  cox_data <- data |> inner_join(clinical.surv, by = "pat_id")
  stopifnot(!(is.na(cox_data$event) |> any()))
  stopifnot(!(is.na(cox_data$time) |> any()))

  results <- map_dfr(score_names, ~ cox.fit(data = cox_data, score_name = .x)) |>
    arrange(p_value)
  results$score_type <- score_type
  results <- results |> mutate(p_adj = p.adjust(p_value, method = "fdr"))

  save.path <- file.path(save_dir, sprintf("figure4_cox_%s.csv", label))
  write_csv(results, save.path)

  pdat <- results |> select(score_name, score_type, hr, hr_lower, hr_upper, p_value, p_adj)
  pdat <- pdat |>
    arrange(p_value) |>
    mutate(
      score_name = factor(score_name, levels = score_name),
      signf = p_adj < 0.05
    )

  p <- ggplot(pdat, aes(y = score_name)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    geom_errorbarh(
      aes(xmin = hr_lower, xmax = hr_upper),
      height = 0.2, linewidth = 0.6, color = "black", alpha = 0.8
    ) +
    geom_point(
      aes(x = hr, fill = signf),
      size = 2.5, shape = 21,
    ) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray"), guide = "none") +
    labs(
      x = "Hazard ratio",
      y = NULL,
      title = paste("Figure 4", label, "-- Cox -", event.name, "~", score_type),
      subtitle = "Points = HR; bars = 95% CI"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 13, face = "bold")
    )
  ggsave(
    filename = file.path(save_dir, sprintf("figure4_cox_%s.png", label)),
    plot = p,
    width = 8,
    height = 6
  )
}

cat("Saved Figure 4d-e panels to", save_dir, "\n")
