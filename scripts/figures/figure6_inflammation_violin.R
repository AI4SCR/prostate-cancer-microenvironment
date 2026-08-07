library(dotenv)
load_dot_env()

library(arrow)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(compositions)
library(rlang)

data_dir <- Sys.getenv("DATA_DIR")
output_figures_dir <- Sys.getenv("OUTPUT_FIGURES_DIR")
stopifnot("DATA_DIR is not set; copy .env.example to .env and fill it in" = nzchar(data_dir))
stopifnot("OUTPUT_FIGURES_DIR is not set; copy .env.example to .env and fill it in" = nzchar(output_figures_dir))

figures_dir <- file.path(output_figures_dir, "figure6")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

compute_label_frequency <- function(data, level, pseudocount = 1) {
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

# Load data
df_clusters <- read_parquet(file.path(data_dir, "niches", "clusters_annotated.parquet"))
df_clusters[["sample_name"]] <- df_clusters[["tma_id"]]

clinical <- read_parquet(file.path(data_dir, "clinical.parquet"))

output_dir <- file.path(output_figures_dir, "figure6")
print(paste("Results will be loaded from / saved to:", output_dir))
print(paste("Figures will be saved to:", figures_dir))



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



# make sure it's a 2-level factor (split violin expects 2 groups)
df$inflammation <-
  factor(df$inflammation,
         levels = c("no", "yes"))

my_cols <- c(
  no   = "#21c8e9",
  yes  = "#e95321"
)


cols_inflammation_up <- c("immune_bloodvessels_CAF1(CD105-)"
                         ,"TLS"
                         ,"Macrophages_Tcells_CAF1(CD105-)")
df_long_full <- df %>%
  pivot_longer(
    cols = all_of(cols_inflammation_up),
    names_to = "niche",
    values_to = "clr_proportion"
  )

my_cols <- c(
  no   = "#21c8e9",
  yes  = "#e95321"
)


p_full <- ggplot(
  df_long_full,
  aes(
    x = inflammation,
    y = clr_proportion,
    fill = inflammation
  )
) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  stat_summary(
    fun.data = "mean_se",
    geom = "pointrange",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = my_cols, name = "Inflammation") +
  facet_wrap(~ niche) +
  labs(
    title = "Inflammatory niches",
    x = "Inflammation",
    y = "CLR-transformed Proportion"
  ) +
  theme_minimal()

# p_full
# ggsave(
#   p_full,
#   filename = file.path(result_dir, "violin_boxplot_inflammation_sep.pdf"),
#   width = 12,
#   height = 18,
#   dpi = 300
# )


pd <- position_dodge(width = 0.9)

p_full <- ggplot(
  df_long_full,
  aes(
    x = niche,
    y = clr_proportion,
    fill = inflammation
  )
) +
  geom_violin(position = pd, alpha = 0.4, trim = FALSE, width = 1.2) +
  geom_boxplot(position = pd, width = 0.25, alpha = 0.6,
               outlier.shape = NA, show.legend = FALSE) +


  scale_fill_manual(values = my_cols, name = "Inflammation") +
  labs(
    title = "Inflammatory niches",
    x = "Niche",
    y = "CLR-transformed Proportion"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# p_full
ggsave(
  p_full,
  filename = file.path(figures_dir, "violin_boxplot_inflammation.pdf"),
  width = 12,
  height = 6,
  dpi = 300
)

cat("Saved Figure 6d panels to", figures_dir, "\n")
