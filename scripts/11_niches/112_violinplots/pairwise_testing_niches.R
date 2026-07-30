library(arrow)
result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies/frequency_boxplots/"

base_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies"
df_props = read_parquet(file.path(base_dir, "stacked_barplots/props_niche_tma_id.parquet"))
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

## make gleason_grp binary
# group 1,2,3 and 4,5
df$gleason_grp_bin <- as.character(df$gleason_grp)
df$gleason_grp_bin[df$gleason_grp_bin %in% c("1", "2", "3")] <- "no"
df$gleason_grp_bin[df$gleason_grp_bin %in% c("4", "5")] <- "yes"
df$gleason_grp_bin <- as.factor(df$gleason_grp_bin)


# target_col <- 'gleason_grp_bin'
# target_col <- "inflammation"
target_col <- "stromogenic_smc_loss_reactive_stroma_present"

## run wilicoxon test for each niche
results <- list()
for (niche in colnames(matrix)){
  group1 <- df %>%
    filter(!!sym(target_col) == 'yes') %>%
    pull(!!sym(niche))
  
  group2 <- df %>%
    filter(!!sym(target_col) == 'no') %>%
    pull(!!sym(niche))
  
  test_result <- wilcox.test(group1, group2)
  
  results[[niche]] <- data.frame(
    niche = niche,
    p_value = test_result$p.value,
    group1_median = median(group1, na.rm = TRUE),
    group2_median = median(group2, na.rm = TRUE)
  )
}
results_df <- do.call(rbind, results)
results_df[['p.adjusted']] <- round(p.adjust(results_df[['p_value']], method = 'BH'), 4)
results_df[['difference_median']] <- results_df[['group1_median']] - results_df[['group2_median']]

# results_df <- results_df %>%
#   arrange(p.adjusted)

## put a star depending on p-value adjusted 
results_df <- results_df %>%
  mutate(significance = case_when(
    p.adjusted < 0.001 ~ '***',
    p.adjusted < 0.01 ~ '**',
    p.adjusted < 0.05 ~ '*',
    TRUE ~ 'ns'
  ))
results_df <- results_df %>%
  mutate(direction = ifelse(difference_median > 0, "pos", "neg"))

save_dir <- paste0(result_dir, "/", target_col, "/")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results_df, file = paste0(save_dir, "pairwise_wilcoxon_results.csv"), row.names = FALSE)



# color mapping for direction
dir_cols <- c(
  pos = "salmon",
  neg = "lightblue"
)

ha_bottom <- HeatmapAnnotation(
  direction = anno_simple(
    results_df$direction,
    col = dir_cols,
    border = FALSE
  ),
  significance = anno_text(
    results_df$significance,
    gp = gpar(col = "black", fontsize = 10),
    just = "center"
  ),
  annotation_height = unit(c(4, 4), "mm"),
  show_annotation_name = FALSE
)

#plot on its own
ht <- Heatmap(
    matrix,
   name = "CLR transformed proportions",
     cluster_rows = TRUE,
     cluster_columns = FALSE,
     show_row_names = FALSE,
     top_annotation = NULL,
     bottom_annotation = ha_bottom,
     column_title = paste0("Niches associated with ", target_col),
   column_title_side = "bottom"
  )
print(ht)




