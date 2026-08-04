library(arrow)
result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches"

base_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/frequencies"
df_props = read_parquet(file.path(base_dir, "stacked_barplots/props_niche_tma_id.parquet"))
# metadata = read_parquet(file.path(base_dir, "metadata_aligned_niche_frequencies.parquet"))
# df_props = read_parquet(file.path(result_dir, "frequencies/stacked_barplots/props_tma.parquet"))
#clusters = read_parquet(file.path(result_dir, "frequencies/stacked_barplots/clustered_metadata.parquet"))
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

clinical.new.path = file.path("/Users/me3312/Downloads/clinical.csv")
clinical = read.csv(clinical.new.path)

tma_ids.valid = intersect(df_props$tma_id, clinical$tma_id)
not_valid = setdiff(clinical$tma_id, tma_ids.valid)
clinical_not_valid = clinical %>%
  filter(tma_id %in% not_valid)

cells <- read_parquet("/Users/me3312/Documents/Paper_PCa/0-paper/0-export/metadata.parquet")
sample_id_valid <- unique(cells$sample_id)
clinical <- clinical %>%
  filter(sample_id %in% sample_id_valid) %>%
  filter(is_tumor == "yes")



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
    d_amico_risk,
    #glandular_atrophy_pin,
    #ln_status,
    tma_id
  ) %>%
  distinct()
 


#df_metadata <- merge(df_metadata, clusters, by = "tma_id", all.x = TRUE)


matrix <- matrix[as.character(df_metadata$tma_id), ]  # reorder matrix rows to match metadata






# hue_list
hue_order_list <- list(
  cause_of_death = c('alive', 'PCa_death', 'non-PCa_death'),
  os_status = c('alive', 'dead'),
  disease_progr = c(0, 1),
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
  ln_status = c("0" = "#2A9D8F", "1" = "#E76F51"),
  os_status = c("alive" = "#06D6A0", "dead" = "#E63946"),
  disease_progr = c("0" = "#457B9D", "1" = "#F4A261")
)


custom_annotation_colors <- list(
  os_status = c(
    alive = "#1f91b4",
    dead  = "#ff0e4a"
  ),
  cause_of_death = c(
    alive        = "#1f91b4",
    PCa_death    = "#ff0e4a",
    `non-PCa_death` = "#ff7f0e"
  ),
  clinical_progr = c(
    "0" = "#1f77b4",
    "1" = "#ffbb78"
  ),
  disease_progr = c(
    "0" = "#1f77b4",
    "1" = "#ffbb78"
  ),
  psa_progr = c(
    "0" = "#1f77b4",
    "1" = "#ffbb78"
  ),
  gs_grp = c(
    "1"   = "#f6d2d2",
    "2"   = "#f1abab",
    "3"   = "#f96b6b",
    "4"   = "#f44336",
    "5"   = "#b71c1c",
    nan   = "#9e9e9e"
  ),
  stromogenic_smc_loss_reactive_stroma_present = c(
    no   = "#9efa70",
    yes  = "#c55797",
    None = "#9e9e9e",
    nan  = "#9e9e9e"
  ),
  inflammation = c(
    no   = "#21c8e9",
    yes  = "#e95321",
    None = "#9e9e9e",
    nan  = "#9e9e9e"
  ),
  # glandular_atrophy_pin = c(
  #   no  = "#1f77b4",
  #   yes = "#ff7f0e",
  #   nan = "#9e9e9e"
  # ),
  glandular_atrophy_pin = c(
    high  = "#b71c1c",
    intermediate = "#ff7f0e",
    nan = "#9e9e9e"
  ),
  gleason_grp = c(
    "1" = "#f6d2d2",
    "2" = "#f1abab",
    "3"= "#f96b6b",
    "4" = "#f44336",
    "5" = "#b71c1c",
    "None"  = "#9e9e9e",
    "nan"   = "#9e9e9e"
  ),
  ln_status = c("0" = "#2A9D8F", "1" = "#E76F51")
)

df_metadata$pat_id <- as.character(df_metadata$pat_id)
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
ha <- rowAnnotation(df = df_metadata, col = annotation_colors, show_legend= c(FALSE, rep(TRUE, length(df_metadata) - 2), FALSE))

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

# mat_clr_transform = as.matrix(compositions::clr(matrix + 1e-7))
# ## denrogram on transformed matrix
# jsd_matrix = as.matrix(dist(mat_clr_transform, method = jsd_distance))
# dend <- as.dendrogram(hclust(as.dist(jsd_matrix), method = "average"))
# dend <- as.dendrogram(hc)
# dend = color_branches(dend, k = 12)
########### -------- create inflammation and stromogenic bottom annotations --------##############
anno_inflammation <- read.csv("/Users/me3312/Documents/Paper_PCa/5-niches/frequencies/frequency_boxplots/inflammation/pairwise_wilcoxon_results.csv")
anno_stromogenic <- read.csv("/Users/me3312/Documents/Paper_PCa/5-niches/frequencies/frequency_boxplots/stromogenic_smc_loss_reactive_stroma_present/pairwise_wilcoxon_results.csv")

anno_inflammation = anno_inflammation %>%
  select(niche, significance, direction)
anno_stromogenic = anno_stromogenic %>%
  select(niche, significance, direction)

df_histo <- merge(anno_inflammation, anno_stromogenic, by = "niche", suffixes = c("_inflammation", "_stromogenic"))

df_histo <- df_histo %>%
  filter(niche %in% colnames(matrix)) %>%
  arrange(match(niche, colnames(matrix)))

# color mapping for direction
dir_cols <- c(
  pos = "salmon",
  neg = "lightblue"
)

ha_bottom <- HeatmapAnnotation(
  inflammation = anno_simple(
    df_histo$direction_inflammation,
    col = dir_cols,
    border = FALSE
  ),
  inflammation_sig = anno_text(
    ifelse(df_histo$significance_inflammation == "ns", "",
           df_histo$significance_inflammation),
    gp = gpar(fontsize = 9, rotation = 90)
  ),
  stromogenic = anno_simple(
    df_histo$direction_stromogenic,
    col = dir_cols,
    border = FALSE
  ),
  stromogenic_sig = anno_text(
    ifelse(df_histo$significance_stromogenic == "ns", "",
           df_histo$significance_stromogenic),
    gp = gpar(fontsize = 9, rotation = 90)
  ),
  annotation_height = unit(c(3, 3, 3, 3), "mm"),
  show_annotation_name = FALSE
)




row_dend = row_dend(ht)
row_dend = color_branches(row_dend, k = 5)

k <- 5
slice_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")  # pick any
# Plot heatmap
p <- ComplexHeatmap::Heatmap(matrix, 
             col = col_fun,
             name = "Proportion", 
             cluster_rows = TRUE, cluster_columns = TRUE,
             clustering_distance_columns = 'pearson',
             clustering_distance_rows = jsd_distance,
             show_row_names = FALSE, show_column_names = TRUE,
             right_annotation = ha,
             top_annotation = niche_anno,
             #bottom_annotation = ha_bottom,
             #column_split = info_niches$meta_niche, # separate visually by meta_niche,
             #row_split = df_metadata$inflammation,
             row_split = 5,
             row_dend_width = unit(3, "cm"),
             #row_dend_gp = lapply(slice_cols, function(cc) gpar(col = cc, lwd = 2)),
             
             heatmap_legend_param = list(title = "Proportion"),
             column_names_rot = 45,
             column_names_gp = gpar(fontsize = 8))
              
# Save heatmap to file
plot_name = "niche_proportion_heatmap_annotated_clustered_final.pdf"
#pdf(paste(base_dir, plot_name, sep = "/"), width = 18, height = 15)
ht <- draw(p)
dev.off()
print(p)



# Extract dendrograms
# 
# row_dend = row_dend(ht)
# hc_row <- as.hclust(row_dend)
# clusters <- cutree(hc_row, k = 5)  # choose number of clusters
# df = data.frame(
#   tma_id = names(clusters),
#   risk_group = as.factor(paste0("cluster_", clusters))
# )


# rownames per cluster
row_orders <- row_order(ht)
clustered_rownames <- lapply(row_orders, \(idx) rownames(matrix)[idx])

# unlist + assign cluster ID based on list index
names(clustered_rownames) <- paste0("cluster_", seq_along(clustered_rownames))
df = stack(clustered_rownames) %>%
  rename(tma_id = values, risk_group = ind)


metadata <- df_metadata
metadata = merge(metadata, df, by = "tma_id", all.x = TRUE)


progression <- clinical %>%
  select(
    tma_id,
    disease_progr_time,
  ) %>%
  distinct()

death <- clinical %>%
  select(
    tma_id,
    os_status,
    last_fu
  ) %>%
  distinct()

death[['overall_survival']] <- ifelse(death[['os_status']] == 'alive', 0, 1)

metadata <- merge(metadata, progression, by = "tma_id", all.x = TRUE)
metadata <- merge(metadata, death, by = "tma_id", all.x = TRUE)

## filter leaves that are less then 10 samples in a cluster
if (!("risk_group" %in% colnames(metadata))) {
  metadata['risk_group'] <- metadata['cluster']
}
metadata_filtered = metadata %>%
  group_by(risk_group) %>%
  filter(n() >= 1) %>%
  ungroup()

metadata_filtered$risk_group <- factor(metadata_filtered$risk_group)

library(survival)
library(survminer)


fit <- survfit(Surv(last_fu, overall_survival) ~ risk_group, data = metadata_filtered)
ggsurvplot(
  fit,
  data = metadata_filtered,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = "Set2",
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Group",
  risk.table.height = 0.25
)
metadata_filtered$disease_progr <- as.numeric(metadata_filtered$disease_progr)
fit2 <- survfit(Surv(disease_progr_time, disease_progr) ~ risk_group, data = metadata_filtered)
ggsurvplot(
  fit2,
  data = metadata_filtered,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = "Set2",
  xlab = "Time",
  ylab = "Progression-free probability",
  legend.title = "Group",
  risk.table.height = 0.25
)



df_samples <- clinical %>%
  select(
    tma_id,
    sample_id
  ) %>%
  distinct() %>%
  inner_join(df, by = "tma_id") %>%
  arrange(risk_group)

#write.csv(df_samples, file.path(result_dir, "frequencies/niche_heatmap_samples_per_cluster_split_12.csv"), row.names = FALSE)


## make table to dataframe and plot as heatmap
tab <- table(metadata$pat_id, metadata$risk_group)
df_tab <- as.data.frame(tab)
df_tab <- df_tab %>%
  rename(
    pat_id = Var1,
    risk_group = Var2,
    count = Freq
  )
df_tab_wide <- df_tab %>%
  pivot_wider(names_from = risk_group, values_from = count, values_fill = 0)
matrix_tab <- df_tab_wide %>%
  column_to_rownames("pat_id") %>%
  as.matrix()

matrix_norm <- matrix_tab / rowSums(matrix_tab)

p_tab <- ComplexHeatmap::Heatmap(matrix_norm, 
                                 col=col_fun,
                             name = "Number of samples", 
                             cluster_rows = TRUE, cluster_columns = FALSE,
                             show_row_names = FALSE, show_column_names = TRUE,
                             heatmap_legend_param = list(title = "Number of samples"),
                             column_names_rot = 45,
                             column_names_gp = gpar(fontsize = 8))

max_per_row <- as.data.frame(apply(matrix_norm, 1, max))
sorted_max_per_row <- max_per_row[order(-max_per_row[,1]), , drop=FALSE]
order <- rownames(sorted_max_per_row)

matrix_norm <- matrix_tab / rowSums(matrix_tab)

matrix_norm <- matrix_norm[order, ]
## make zeors to na
#matrix_norm[matrix_norm == 0] <- NA
col_fun <- colorRamp2(c(0, 1), c("white", "darkred"))

df_patient <- df_metadata %>%
  select(
    pat_id,
    disease_progr,
    gs_grp,
    os_status,
  ) %>%
  distinct()
#order by order
df_patient <- df_patient[match(rownames(matrix_norm), df_patient$pat_id), ]

right_anno <- rowAnnotation(df = df_patient, col = list(
  disease_progr = c("0" = "#457B9D", "1" = "#E63946"),
  gs_grp = c("1" = "yellow", "2" = "#E9C46A", "3" = "#F4A261", "4" = "#E76F51", "5" = "red", "nan" = "#999999"),
  os_status = c("alive" = "#06D6A0", "dead" = "#E63946")
), show_legend= c(FALSE, TRUE, TRUE, TRUE))

p_tab <- ComplexHeatmap::Heatmap(matrix_norm, 
                                 col=col_fun,
                                 #na_col = "white",
                                 name = "Number of samples", 
                                 cluster_rows = FALSE, cluster_columns = TRUE,
                                 show_row_names = FALSE, show_column_names = TRUE,
                                 heatmap_legend_param = list(title = "Number of samples"),
                                 rect_gp = gpar(col = "black", lwd = 0.1),
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize = 8),
                                 right_annotation = right_anno)
print(p_tab)
