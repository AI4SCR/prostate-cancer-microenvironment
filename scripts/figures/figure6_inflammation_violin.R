# Reproduce Figure 6d: niche composition (CLR-transformed) vs inflammation.
#
# 1:1 port of the old repo's
# 000_paper/11_niches/112_violinplots/inflammation_vis.R (paths changed,
# unused `ggpubr`/`entropy` imports dropped; see figure_script_mapping.md).
# clusters_annotated_v2.parquet has no reproducing script in this repo (see
# figure5_niche_annotation.py's docstring), so it's read from
# LEGACY_DATA_DIR's precomputed copy; clinical.parquet comes from EXPORT_DIR.
#
# The original script's first two panels use `introdataviz::geom_split_violin`
# (a GitHub-only package, not on CRAN) for a split-by-inflammation violin.
# Neither `introdataviz` nor the CRAN alternative `gghalves` install in this
# environment (`gghalves` isn't built for this R version; no `remotes`/
# `devtools` here to install from GitHub). Substituted with plain
# `geom_violin()` + `position_dodge()` -- the SAME visual approach the
# original script already uses for its other two panels on the same data
# (see its "sep"/dodge variants), just applied consistently to all three
# instead of split-violins for the first two. Same data, same statistics,
# side-by-side instead of split violins.
#
# Writes to $EXPORT_DIR/figures/figure6/.

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
stopifnot("EXPORT_DIR is not set; copy .env.example to .env and fill it in" = nzchar(export_dir))
stopifnot("LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(legacy_dir))

save_dir <- file.path(export_dir, "figures", "figure6")
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
df_long$inflammation <- factor(df_long$inflammation, levels = c("no", "yes"))

my_cols <- c(no = "#21c8e9", yes = "#e95321")
pd <- position_dodge(width = 0.9)

p_up <- ggplot(df_long, aes(x = niche, y = clr_proportion, fill = inflammation)) +
  geom_violin(position = pd, alpha = 0.4, trim = FALSE) +
  geom_boxplot(position = pd, width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  facet_wrap(~direction, scales = "free_y") +
  labs(title = "inflammation niches upregulated in reactive stroma", x = "Niche", y = "CLR-transformed Proportion") +
  theme_minimal()
ggsave(file.path(save_dir, "figure6d_inflammation_updown_niches.pdf"), p_up, width = 12, height = 6, dpi = 300)

df_long_full <- df %>%
  pivot_longer(cols = all_of(colnames(df_clr)[-1]), names_to = "niche", values_to = "clr_proportion")

p_split <- ggplot(df_long_full, aes(x = niche, y = clr_proportion, fill = inflammation)) +
  geom_violin(position = pd, alpha = 0.4, trim = FALSE) +
  geom_boxplot(position = pd, width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  labs(title = "All niches", x = "Niche", y = "CLR-transformed Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(save_dir, "figure6d_inflammation_all_niches_split.pdf"), p_split, width = 12, height = 6, dpi = 300)

p_sep <- ggplot(df_long_full, aes(x = inflammation, y = clr_proportion, fill = inflammation)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = FALSE) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  facet_wrap(~niche) +
  labs(title = "All niches", x = "Reactive stroma", y = "CLR-transformed Proportion") +
  theme_minimal()
ggsave(file.path(save_dir, "figure6d_inflammation_all_niches_sep.pdf"), p_sep, width = 12, height = 18, dpi = 300)

cat("Saved Figure 6d panels to", save_dir, "\n")
