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



#df_metadata <- merge(df_metadata, clusters, by = "tma_id", all.x = TRUE)


matrix <- matrix[as.character(df_metadata$tma_id), ]  # reorder matrix rows to match metadata

info_niches = read.csv(file.path(result_dir, "/annotation/niche_annotations_v2.csv"))

info_niches = info_niches %>%
  select(-cluster) %>%
  distinct() %>%
  #filter(niche %in% rownames(matrix)) %>%
  filter(niche != 'unassigned')


# reorder matrix cols by niche
matrix = matrix[, info_niches$niche]

### 
# calculate correlation between all columns of matrix
corr_matrix <- cor(matrix, use = "pairwise.complete.obs", method = "spearman")
## plot as heatmap with numbers in it

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#2166ac", "white", "#b2182b")
)

h <- Heatmap(
  corr_matrix,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(
      sprintf("%.2f", corr_matrix[i, j]),
      x, y,
      gp = grid::gpar(fontsize = 8)
    )
  }
)
plot_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/visualization/composition"
plot_name = file.path(plot_dir, "niche_correlation_spearman_heatmap_numbers.pdf")

pdf(plot_name, width = 18, height = 14)
draw(h, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()



