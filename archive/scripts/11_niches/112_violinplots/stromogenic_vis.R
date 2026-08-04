library(arrow)
result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies/frequency_boxplots/"

base_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies"
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(entropy)
library(compositions)

clinical.path = file.path('/Users/me3312/Documents/Paper_PCa/0-paper/0-export/clinical.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()




compute_label_frequency <- function(data, level, pseudocount = 0) {
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
df_freqs <- compute_label_frequency(df_clusters, level = "niche", pseudocount = 1)

## rename sample_name to tma_id, select niche, proportion and make to wide format
df_wide <- df_freqs %>%
  select(tma_id = sample_name, niche, proportion) %>%
  pivot_wider(names_from = niche, values_from = proportion, values_fill = 0)

df_props <- df_wide

tma_ids.valid = intersect(df_props$tma_id, clinical$tma_id)


matrix = df_props %>%
  column_to_rownames("tma_id") %>%
  as.matrix()

clinical <- clinical %>%
  filter(tma_id %in% rownames(matrix)) %>%
  arrange(tma_id)

# clusters <- clusters %>%
#   select(tma_id, cluster) %>%
#   distinct()

df_metadata <- clinical %>%
  select(
    pat_id,
    os_status,
    disease_progr,
    #psa_progr,
    #recurrence,
    gs_grp,
    gleason_grp,
    inflammation, 
    stromogenic_smc_loss_reactive_stroma_present,
    #glandular_atrophy_pin,
    #ln_status,
    tma_id
  ) %>%
  distinct()

df2 <- column_to_rownames(df_props, var = "tma_id")

# 2. CLR transform (row-wise)
df_clr <- as.data.frame(clr(df2))


# 3. Put ID back as column
df_clr <- rownames_to_column(df_clr, var = "tma_id")

df <- df_clr %>%
  inner_join(df_metadata, by = "tma_id")

cols_stromogenic_up <- c("tumor_CAF1(CD105High)"
                         ,"tumor_CAF1_lymphocytes"
                         ,"luminal_CAF1(CD105High)")
  
cols_stromogenic_down <- c("CAF2(AR-)_enriched",
                           "CAF2s_enriched",
                           "CAFs_lymphocytes")

df_long <- df %>%
  pivot_longer(
    cols = all_of(c(cols_stromogenic_up, cols_stromogenic_down)),
    names_to = "niche",
    values_to = "clr_proportion"
  )
df_long$direction <- ifelse(df_long$niche %in% cols_stromogenic_up, "up", "down")


library(ggplot2)
library(introdataviz)

# make sure it's a 2-level factor (split violin expects 2 groups)
df_long$stromogenic_smc_loss_reactive_stroma_present <-
  factor(df_long$stromogenic_smc_loss_reactive_stroma_present,
         levels = c("no", "yes"))

my_cols <- c(
  no   = "#9efa70",
  yes  = "#c55797",
)

p_up <- ggplot(
  df_long,
  aes(
    x = niche,
    y = clr_proportion,
    fill = stromogenic_smc_loss_reactive_stroma_present
  )
) +
  introdataviz::geom_split_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    show.legend = FALSE,
    position = position_dodge(0.175)
  ) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  facet_wrap(~ direction, scales = "free_y") +
  labs(
    title = "Stromogenic niches upregulated in reactive stroma",
    x = "Niche",
    y = "CLR-transformed Proportion"
  ) +
  theme_minimal()

p_up

result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies/frequency_boxplots/violin_stromogenic/"
dir.create(result_dir, showWarnings = FALSE)

df_long_full <- df %>%
  pivot_longer(
    cols = all_of(colnames(df_clr)[-1]),
    names_to = "niche",
    values_to = "clr_proportion"
  )

my_cols <- c(
  no   = "#9efa70",
  yes  = "#c55797"
)
p_full <- ggplot(
  df_long_full,
  aes(
    x = niche,
    y = clr_proportion,
    fill = stromogenic_smc_loss_reactive_stroma_present
  )
) +
  introdataviz::geom_split_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    show.legend = FALSE,
    position = position_dodge(0.175)
  ) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  labs(
    title = "All niches",
    x = "Niche",
    y = "CLR-transformed Proportion"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_full
ggsave(
  p_full,
  filename = file.path(result_dir, "violin_boxplot_stromogenic_split.pdf"),
  width = 12,
  height = 6,
  dpi = 300
)

p_full <- ggplot(
  df_long_full,
  aes(
    x = stromogenic_smc_loss_reactive_stroma_present,
    y = clr_proportion,
    fill = stromogenic_smc_loss_reactive_stroma_present
  )
) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  facet_wrap(~ niche) +
  labs(
    title = "All niches",
    x = "Reactive stroma",
    y = "CLR-transformed Proportion"
  ) +
  theme_minimal()

p_full
ggsave(
  p_full,
  filename = file.path(result_dir, "violin_boxplot_stromogenic_sep.pdf"),
  width = 12,
  height = 18,
  dpi = 300
)


pd <- position_dodge(width = 0.9)

p_full <- ggplot(
  df_long_full,
  aes(
    x = niche,
    y = clr_proportion,
    fill = stromogenic_smc_loss_reactive_stroma_present
  )
) +
  geom_violin(position = pd, alpha = 0.4, trim = FALSE, width = 1.2) +
  geom_boxplot(position = pd, width = 0.25, alpha = 0.6,
               outlier.shape = NA, show.legend = FALSE) +
  

  scale_fill_manual(values = my_cols, name = "Reactive stroma") +
  labs(
    title = "All niches",
    x = "Niche",
    y = "CLR-transformed Proportion"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_full
ggsave(
  p_full,
  filename = file.path(result_dir, "violin_boxplot_stromogenic.pdf"),
  width = 12,
  height = 6,
  dpi = 300
)

