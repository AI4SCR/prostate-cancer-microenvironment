# Reproduce Figure 7a: niche composition (CLR-transformed) vs stromogenic
# status.
#
# 1:1 port of the old repo's
# 000_paper/11_niches/112_violinplots/stromogenic_vis.R (paths changed,
# unused `ComplexHeatmap`/`circlize`/`ggpubr`/`entropy` imports dropped).
# clusters_annotated_v2.parquet has no reproducing script in this repo, so
# it's read from LEGACY_DATA_DIR's precomputed copy; clinical.parquet comes
# from EXPORT_DIR. Does not touch figure6_inflammation_violin.R (its
# sibling script, same structure, different clinical variable) or the
# already-validated figure6_km_niche6.R family -- new, dedicated script.
#
# Same `introdataviz::geom_split_violin` substitution as
# figure6_inflammation_violin.R -- neither `introdataviz` nor `gghalves`
# install in this environment; substituted with plain `geom_violin()` +
# `position_dodge()`, the same visual approach legacy's own other panels
# already use on the same data.
#
# Disclosed fix: legacy's `my_cols <- c(no = "#9efa70", yes = "#c55797",)`
# has a trailing comma inside `c(...)`, which is not valid R syntax --
# confirmed this makes the legacy script unparseable past this line as
# literally written. Fixed by removing the trailing comma (a syntax-only
# fix, not a logic change).
#
# Legacy's `p_up` panel (stromogenic-specific up/down niches, faceted) is
# built and printed but never saved via ggsave -- unlike its
# `figure6_inflammation_violin.R` sibling, where the equivalent panel IS
# saved. Matched verbatim (not forced to save) since legacy already
# produces real, meaningful output here (the 3 "all niches" panels below);
# this isn't the "computed but the whole panel is otherwise never produced"
# case that warranted enabling a save elsewhere in this project (e.g.
# Figure 6e, Supplementary Fig 5b/6a).
#
# Writes to $OUTPUT_FIGURES_DIR/figure7/.

library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(compositions)
library(rlang)

export_dir <- Sys.getenv("EXPORT_DIR")
legacy_dir <- Sys.getenv("LEGACY_DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

save_dir <- file.path(output_figures_dir, "figure7")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- read_parquet(file.path(export_dir, "clinical.parquet"))
num.patients <- clinical$pat_id |> n_distinct()

compute_label_frequency <- function(data, level, pseudocount = 0) {
  level_sym <- rlang::sym(level)
  data_summary <- data %>%
    dplyr::group_by(sample_name, !!level_sym) %>%
    dplyr::summarise(count = dplyr::n() + pseudocount, .groups = "drop") %>%
    tidyr::complete(sample_name, !!level_sym, fill = list(count = pseudocount))
  data_summary %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(proportion = count / sum(count)) %>%
    dplyr::ungroup()
}

df_clusters <- read_parquet(file.path(legacy_dir, "5-niches", "annotation", "clusters_annotated_v2.parquet"))
df_clusters[["sample_name"]] <- df_clusters[["tma_id"]]

df_freqs <- compute_label_frequency(df_clusters, level = "niche", pseudocount = 1)
df_props <- df_freqs %>%
  select(tma_id = sample_name, niche, proportion) %>%
  pivot_wider(names_from = niche, values_from = proportion, values_fill = 0)

matrix <- df_props %>% column_to_rownames("tma_id") %>% as.matrix()
clinical <- clinical %>% filter(tma_id %in% rownames(matrix)) %>% arrange(tma_id)

df_metadata <- clinical %>%
  select(pat_id, os_status, disease_progr, gs_grp, gleason_grp, inflammation, stromogenic_smc_loss_reactive_stroma_present, tma_id) %>%
  distinct()

df_clr <- as.data.frame(clr(column_to_rownames(df_props, var = "tma_id")))
df_clr <- rownames_to_column(df_clr, var = "tma_id")
df <- df_clr %>% inner_join(df_metadata, by = "tma_id")

cols_stromogenic_up <- c("tumor_CAF1(CD105High)", "tumor_CAF1_lymphocytes", "luminal_CAF1(CD105High)")
cols_stromogenic_down <- c("CAF2(AR-)_enriched", "CAF2s_enriched", "CAFs_lymphocytes")

df_long <- df %>%
  pivot_longer(cols = all_of(c(cols_stromogenic_up, cols_stromogenic_down)), names_to = "niche", values_to = "clr_proportion")
df_long$direction <- ifelse(df_long$niche %in% cols_stromogenic_up, "up", "down")
df_long$stromogenic_smc_loss_reactive_stroma_present <- factor(df_long$stromogenic_smc_loss_reactive_stroma_present, levels = c("no", "yes"))

my_cols <- c(no = "#9efa70", yes = "#c55797")
pd <- position_dodge(width = 0.9)

p_up <- ggplot(df_long, aes(x = niche, y = clr_proportion, fill = stromogenic_smc_loss_reactive_stroma_present)) +
  geom_violin(position = pd, alpha = 0.4, trim = FALSE) +
  geom_boxplot(position = pd, width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = FALSE, position = position_dodge(0.175)) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  facet_wrap(~direction, scales = "free_y") +
  labs(title = "Stromogenic niches upregulated in reactive stroma", x = "Niche", y = "CLR-transformed Proportion") +
  theme_minimal()
# p_up computed but not saved, matching legacy verbatim (see docstring)

df_long_full <- df %>%
  pivot_longer(cols = all_of(colnames(df_clr)[-1]), names_to = "niche", values_to = "clr_proportion")

p_split <- ggplot(df_long_full, aes(x = niche, y = clr_proportion, fill = stromogenic_smc_loss_reactive_stroma_present)) +
  geom_violin(position = pd, alpha = 0.4, trim = FALSE) +
  geom_boxplot(position = pd, width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = FALSE, position = position_dodge(0.175)) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  labs(title = "All niches", x = "Niche", y = "CLR-transformed Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(save_dir, "figure7a_stromogenic_all_niches_split.pdf"), p_split, width = 12, height = 6, dpi = 300)

p_sep <- ggplot(df_long_full, aes(x = stromogenic_smc_loss_reactive_stroma_present, y = clr_proportion, fill = stromogenic_smc_loss_reactive_stroma_present)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = FALSE) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  facet_wrap(~niche) +
  labs(title = "All niches", x = "Reactive stroma", y = "CLR-transformed Proportion") +
  theme_minimal()
ggsave(file.path(save_dir, "figure7a_stromogenic_all_niches_sep.pdf"), p_sep, width = 12, height = 18, dpi = 300)

p_full <- ggplot(df_long_full, aes(x = niche, y = clr_proportion, fill = stromogenic_smc_loss_reactive_stroma_present)) +
  geom_violin(position = pd, alpha = 0.4, trim = FALSE, width = 1.2) +
  geom_boxplot(position = pd, width = 0.25, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  labs(title = "All niches", x = "Niche", y = "CLR-transformed Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(save_dir, "figure7a_stromogenic_all_niches.pdf"), p_full, width = 12, height = 6, dpi = 300)

cat("Saved Figure 7a panels to", save_dir, "\n")
