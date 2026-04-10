library(arrow)
result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches"

base_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies"
df_props = read_parquet(file.path(base_dir, "stacked_barplots/props_niche_pat_id.parquet"))
# metadata = read_parquet(file.path(base_dir, "metadata_aligned_niche_frequencies.parquet"))
# df_props = read_parquet(file.path(result_dir, "frequencies/stacked_barplots/props_tma.parquet"))
# clusters = read_parquet(file.path(result_dir, "frequencies/stacked_barplots/clustered_metadata.parquet"))
# ln_status = read_parquet(file.path(result_dir, "frequencies/ln_status.parquet"))
# clinical = read_parquet(file.path(result_dir, "frequencies/clinical.parquet"))
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(entropy)

clinical.path = file.path('/Users/me3312/Documents/Paper_PCa/0-paper/0-export/clinical.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

##save clinical to desktop
 #write.csv(clinical, file = '/Users/me3312/Desktop/clinical.csv', row.names = FALSE)

tma_ids.valid = intersect(df_props$pat_id, clinical$pat_id)

## order metadata and df_probs by tma_id

# metadata <- merge(metadata, ln_status, by = "tma_id", all.x = TRUE)
# metadata = metadata %>%
#   arrange(tma_id)
# df_props = df_props %>%
#   arrange(tma_id)
# metadata$tma_id == df_props$tma_id  # check order matches

matrix = df_props %>%
  column_to_rownames("pat_id") %>%
  as.matrix()

clinical <- clinical %>%
  filter(pat_id %in% rownames(matrix)) %>%
  arrange(pat_id)


df_metadata <- clinical %>%
  select(
    cause_of_death,
    clinical_progr,
    psa_progr,
    recurrence,
    gs_grp,
    pat_id
  ) %>%
  distinct()



matrix <- matrix[as.character(df_metadata$pat_id), ]  # reorder matrix rows to match metadata






# hue_list
hue_order_list <- list(
  cause_of_death = c('alive', 'PCa_death', 'non-PCa_death'),
  clinical_progr = c(0, 1),
  psa_progr = c(0, 1),
  recurrence = c('no_recurrence', 'recurrence'),
  gs_grp = c('1', '2', '3', '4', '5', 'nan'),
  gleason_grp = c('1', '2', '3', '4', '5'), 
  inflammation = c('yes', 'no', 'nan'),
  stromogenic_smc_loss_reactive_stroma_present = c('yes', 'no', 'nan'),
  glandular_atrophy_pin = c('yes', 'no', 'nan'),
  #cluster = seq(0, length(unique(df_metadata$cluster)), 1),
  ln_status = c(0, 1)
)



# Define annotation colors (auto-detect categorical columns)
for (col in names(hue_order_list)) {
  if (col %in% colnames(df_metadata)) {
    df_metadata[[col]] <- factor(df_metadata[[col]], levels = hue_order_list[[col]])
  }
}


# Define custom colors for each metadata column
custom_annotation_colors <- list(
  clinical_progr = c("0" = "#457B9D", "1" = "#E63946"),
  cause_of_death = c("alive" = "#06D6A0", "PCa_death" = "#EF476F", "non-PCa_death" = "#457B9D"),
  recurrence = c("no_recurrence" = "#1D3557", "recurrence" = "#F4A261"),
  gleason_grp = c("1" = "yellow", "2" = "#E9C46A", "3" = "#F4A261", "4" = "#E76F51", "5" = "red"),
  inflammation = c("yes" = "#EF476F", "no" = "#457B9D", "nan" = "#999999"),
  stromogenic_smc_loss_reactive_stroma_present = c("yes" = "#E76F51", "no" = "#2A9D8F", "nan" = "#999999"),
  glandular_atrophy_pin = c("yes" = "#F4A261", "no" = "#264653", "nan" = "#999999"),
  ln_status = c("0" = "#2A9D8F", "1" = "#E76F51")
)


annotation_colors <- list()
for (col in colnames(df_metadata)) {
  if (col %in% names(custom_annotation_colors)) {
    annotation_colors[[col]] <- custom_annotation_colors[[col]]
  } else if (col %in% names(hue_order_list)) {
    # Generate consistent colors based on predefined order
    unique_vals <- hue_order_list[[col]]
    annotation_colors[[col]] <- structure(
      circlize::rand_color(length(unique_vals)), 
      names = unique_vals
    )
  }
}

# Create HeatmapAnnotation
ha <- rowAnnotation(df = df_metadata, col = annotation_colors, show_legend= c(rep(TRUE, 5), FALSE))

col_fun <- colorRamp2(c(0, 1), c("white", "darkgreen"))

jsd_distance <- function(x, y){
  entropy_x <- entropy(x)
  entropy_y <- entropy(y)
  m <- (x + y) / 2
  
  ### JSD
  jsd <- entropy(m) - (entropy_x + entropy_y) / 2
  return(jsd)
  
}





# --- define annotation colors ---

info_niches = read.csv(file.path(result_dir, "/annotation/niche_annotations_v2.csv"))

info_niches = info_niches %>%
  select(-cluster) %>%
  distinct() %>%
  #filter(niche %in% rownames(matrix)) %>%
  filter(niche != 'unassigned')


# reorder matrix cols by niche
matrix = matrix[, info_niches$niche]


# build named color mappings (unique)
niche_colors <- setNames(unique(info_niches$niche_color), unique(info_niches$niche))
meta_niche_colors <- setNames(unique(info_niches$meta_niche_color), unique(info_niches$meta_niche))

# build annotation as a data frame
niche_anno <- HeatmapAnnotation(
  df = data.frame(
    meta_niche = factor(info_niches$meta_niche, levels = names(meta_niche_colors)),
    niche = factor(info_niches$niche, levels = names(niche_colors))
  ),
  col = list(
    meta_niche = meta_niche_colors,
    niche = niche_colors
  ),
  show_legend = c(TRUE, FALSE)
)
# 
# mat_clr_transform = as.matrix(compositions::clr(matrix + 1e-7))
# ## denrogram on transformed matrix
# jsd_matrix <- as.matrix(distance(matrix, method = "jensen-shannon"))
# dend <- as.dendrogram(hclust(as.dist(jsd_matrix), method = "average"))
# dend = color_branches(dend, k = 12)
# plot(dend)


# Plot heatmap
p <- ComplexHeatmap::Heatmap(matrix, 
             col = col_fun,
             name = "Proportion", 
             cluster_rows = TRUE, cluster_columns = TRUE,
             clustering_distance_columns = 'spearman',
             clustering_distance_rows = jsd_distance,
             show_row_names = FALSE, show_column_names = TRUE,
             right_annotation = ha,
             top_annotation = niche_anno,
             column_split = info_niches$meta_niche, # separate visually by meta_niche,
             #row_split = df_metadata$inflammation,
             row_split = 5,
             heatmap_legend_param = list(title = "Proportion"),
             column_names_rot = 45,
             column_names_gp = gpar(fontsize = 8))
# Save heatmap to file
plot_name = "niche_proportion_heatmap_annotated_patient.pdf"
#pdf(paste(base_dir, plot_name, sep = "/"), width = 28, height = 20)
ht <- draw(p)
dev.off()
print(p)



# Extract dendrograms
# 
# row_dend = row_dend(ht)
# hc_row <- as.hclust(row_dend)
# clusters <- cutree(hc_row, k = 5)  # choose number of clusters
# df = data.frame(
#   pat_id = names(clusters),
#   risk_group = as.factor(paste0("cluster_", clusters))
# )


# rownames per cluster
row_orders <- row_order(ht)
clustered_rownames <- lapply(row_orders, \(idx) rownames(matrix)[idx])

# unlist + assign cluster ID based on list index
names(clustered_rownames) <- paste0("cluster_", seq_along(clustered_rownames))
df = stack(clustered_rownames) %>%
  rename(pat_id = values, risk_group = ind)


metadata <- df_metadata
metadata = merge(metadata, df, by = "pat_id", all.x = TRUE)


progression <- clinical %>%
  select(
    pat_id,
    clinical_progr_time,
  ) %>%
  distinct()

death <- clinical %>%
  select(
    pat_id,
    os_status,
    last_fu
  ) %>%
  distinct()

death[['overall_survival']] <- ifelse(death[['os_status']] == 'alive', 0, 1)

metadata <- merge(metadata, progression, by = "pat_id", all.x = TRUE)
metadata <- merge(metadata, death, by = "pat_id", all.x = TRUE)

## filter leaves that are less then 10 samples in a cluster
if (!("risk_group" %in% colnames(metadata))) {
  metadata['risk_group'] <- metadata['cluster']
}
metadata_filtered = metadata %>%
  group_by(risk_group) %>%
  filter(n() >= 10) %>%
  ungroup()

metadata_filtered$risk_group <- factor(metadata_filtered$risk_group)
metadata_filtered$clinical_progr <- as.numeric(metadata_filtered$clinical_progr)
  
library(survival)
library(survminer)
fit <- survfit(Surv(clinical_progr_time, clinical_progr) ~ risk_group, data = metadata_filtered)
ggsurvplot(
  fit,
  data = metadata_filtered,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = "Set3",
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Group",
  risk.table.height = 0.25
)


fit <- survfit(Surv(last_fu, overall_survival) ~ risk_group, data = metadata_filtered)
ggsurvplot(
  fit,
  data = metadata_filtered,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = "Set3",
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Group",
  risk.table.height = 0.25
)

df_patient <- clinical %>%
  select(
    pat_id,
    os_status,
    cause_of_death,
    clinical_progr,
    clinical_progr_time,
    last_fu
  ) %>%
  distinct()
