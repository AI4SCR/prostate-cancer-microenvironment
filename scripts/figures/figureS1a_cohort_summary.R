library(dotenv)
load_dot_env()

library(arrow)
library(tidyverse)
library(rlang)

export_dir <- Sys.getenv("EXPORT_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))

outputs_dir <- file.path(output_figures_dir, "figureS1")

get_stacked_barplot <- function(data, var_name) {
  data[[var_name]] <- factor(data[[var_name]])

  counts <- data %>%
    count(.data[[var_name]]) %>%
    mutate(label = paste0(.data[[var_name]], " (n=", n, ")"))

  labels_map <- setNames(counts$label, counts[[var_name]])

  g <- ggplot(data, aes(x = 1, fill = .data[[var_name]])) +
    geom_bar(position = "stack") +
    scale_fill_discrete(labels = labels_map) +
    labs(x = NULL, y = "Number of observations", fill = var_name) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.title = element_blank()
    )

  return(g)
}

get_violin_plot <- function(data, var_name) {
  g <- ggplot(data, aes(x = 1, y = .data[[var_name]])) +
    geom_violin(width = 1, trim = FALSE, alpha = 0.4, fill = "grey80", color = "grey50") +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(width = 0.125, size = 1, alpha = 0.6, color = "black") +
    labs(x = NULL, y = var_name) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.title = element_blank()
    )
  g
  return(g)
}

set.seed(1)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
clinical <- clinical[order(clinical$sample_id), ]

filter_ <- duplicated(clinical$tma_id)
tmas <- clinical[!filter_, ]

filter_ <- duplicated(clinical$pat_id)
pats <- clinical[!filter_, ]

# %% TMA-level variables
tma.cat_names <- c(
  "gs_pat_1",
  "gs_pat_2",
  "gleason_grp",
  "stromogenic_smc_loss_reactive_stroma_present",
  "non_stromogenic_smc_abundant",
  "inflammation",
  "glandular_atrophy_pin",
  "cribriform",
  "is_tumor",
  "gleason_score",
  "gleason_score_sum",
  "gleason_grp"
)

# %% PATIENT-level variables
pat.cat_names <- c(
  "cause_of_death",
  "os_status",
  "psa_progr",
  "clinical_progr",
  "disease_progr",
  "recurrence_loc",
  "cgrading_biopsy",
  "cgs_score",
  "gs_grp",
  "ct_stage",
  "pt_stage",
  "pgs_score",
  "ln_status",
  "surgical_margin_status",
  "adj_adt",
  "adj_radio",
  "recurrence",
  "d_amico_risk",
  "os_status",
  "pt_stage"
)

pat.con_names <- c(
  "age_at_surgery",
  "psa_at_surgery",
  "last_fu",
  "psa_progr_time",
  "clinical_progr_time",
  "disease_progr_time"
)

# %% CATEGORICAL
for (var_name in tma.cat_names) {
  g <- get_stacked_barplot(tmas, var_name)

  fname <- paste0(var_name, ".pdf")
  save_path <- file.path(outputs_dir, "tma-level", fname)
  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, plot = g, width = 3, height = 4)
}

for (var_name in pat.cat_names) {
  g <- get_stacked_barplot(pats, var_name)

  fname <- paste0(var_name, ".pdf")
  save_path <- file.path(outputs_dir, "patient-level", fname)
  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, plot = g, width = 3, height = 4)
}

# %% CONTINUOUS
for (var_name in pat.con_names) {
  # only plot patients with events
  if (var_name == "psa_progr_time") {
    data <- pats %>% filter(psa_progr == 1)
  } else if (var_name == "clinical_progr_time") {
    data <- pats %>% filter(clinical_progr == 1)
  } else if (var_name == "disease_progr_time") {
    data <- pats %>% filter(disease_progr == 1)
  } else {
    data <- pats
  }

  g <- get_violin_plot(data = data, var_name = var_name)

  fname <- paste0(var_name, ".pdf")
  save_path <- file.path(outputs_dir, "patient-level", fname)
  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, plot = g, width = 3, height = 4)
}

cat("Saved Supplementary Figure 1a panels to", outputs_dir, "\n")
