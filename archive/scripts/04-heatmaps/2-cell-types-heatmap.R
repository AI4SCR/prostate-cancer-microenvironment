# install.packages("arrow")
# install.packages("circlize")
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("ComplexHeatmap")

library(arrow)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(yaml)
library(viridis)
library(stringr)

save_dir = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/outputs/3-heatmaps')
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
data_path = file.path(
  "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/intensity_normalized.parquet"
)
metadata_path = file.path(
  "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/metadata.parquet"
)
colormap_path = file.path("/Users/adrianomartinelli/projects/PCa/colormaps.yaml")
index_names = c('sample_id', 'object_id')

main_group = 'immune'
subset_size = 2000
cluster_rows = FALSE
exclude_channels = c('fap', 'dna1', 'dna2', 'icsk1', 'icsk2', 'icsk3')

heatmap = function(data_path,
                   metadata_path,
                   colormap_path,
                   index_names = c('sample_id', 'object_id'),
                   cluster_rows = FALSE,
                   subset_size = 2000,
                   aggregate_by = NULL) {
  # -----------------------------------------------------------------------------
  # LOAD DATA
  data <- read_parquet(data_path)
  data = data |> select(-all_of(exclude_channels)) 
  
  meta <- read_parquet(metadata_path)
  colormaps <- yaml::read_yaml(colormap_path)
  
  # -----------------------------------------------------------------------------
  # SUBSET (only if not aggregating, otherwise keep all)
  if (is.null(aggregate_by)) {
    set.seed(0)
    subset = sample(nrow(data), subset_size)
    data = data[subset, ]
    meta = meta[subset, ]
  }
  
  # -----------------------------------------------------------------------------
  # SORT
  ord <- order(meta$label)
  
  data <- data[ord, ]
  meta <- meta[ord, ]
  
  # -----------------------------------------------------------------------------
  # INDEX
  data.index = data[, index_names]
  data = data[, !(names(data) %in% index_names)]
  
  meta.index = meta[, index_names]
  meta = meta[, !(names(meta) %in% index_names)]
  
  stopifnot(all(meta.index == data.index))
  index = do.call(paste0, meta.index[index_names])
  
  # order = match(df1, rownames(df2))
  
  mat <- as.matrix(data)
  rownames(mat) <- index
  
  # -----------------------------------------------------------------------------
  # AGGREGATION (mean per group)
  if (!is.null(aggregate_by)) {
    group <- meta[[aggregate_by]]
    mat <- rowsum(mat, group) / as.vector(table(group))
    meta <- meta[!duplicated(group), , drop = FALSE]
    stopifnot(all(meta[[aggregate_by]] == rownames(mat)))
  }
  
  # Define colour mappings for annotations
  label_colors <- unlist(colormaps$label)
  main_group_colors <- unlist(colormaps$main_group)
  patient_colors <- unlist(colormaps$pat_id)
  
  # -----------------------------------------------------------------------------
  # Create row annotations
  col_anno <- HeatmapAnnotation(
    label      = meta$label,
    main_group = meta$main_group,
    patient = meta$pat_id,
    col = list(
      label      = label_colors,
      main_group = main_group_colors,
      patient    = patient_colors
    ),
    na_col = "#F0F0F0",
    show_legend = c(
      label = FALSE,
      main_group = TRUE,
      patient = FALSE
    )
  )
  
  # -----------------------------------------------------------------------------
  # Draw the heatmap
  heatmap_obj <- Heatmap(
    t(mat),   # transpose = rotate 90 degrees
    name = "Protein intensity",
    col = inferno(256),
    cluster_rows = TRUE,   # swap
    cluster_columns = cluster_rows,   # swap
    show_row_names = TRUE,            # was columns
    show_column_names = FALSE,        # was rows
    column_split = meta$main_group,   # was row_split
    column_title = NULL,              # was row_title
    top_annotation = col_anno,        # was right_annotation
    heatmap_legend_param = list(title = "Intensity")
  )
  
  draw(heatmap_obj,
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
}

heatmap.agg = function(data_path,
                   metadata_path,
                   colormap_path,
                   index_names = c('sample_id', 'object_id'),
                   cluster_rows = FALSE,
                   subset_size = 2000,
                   aggregate_by = NULL,
                   column_split = TRUE) {
  # -----------------------------------------------------------------------------
  # LOAD DATA
  data <- read_parquet(data_path)
  data = data |> select(-all_of(exclude_channels)) 
  
  meta <- read_parquet(metadata_path)
  colormaps <- yaml::read_yaml(colormap_path)
  
  # -----------------------------------------------------------------------------
  # SUBSET (only if not aggregating, otherwise keep all)
  if (is.null(aggregate_by)) {
    set.seed(0)
    subset = sample(nrow(data), subset_size)
    data = data[subset, ]
    meta = meta[subset, ]
  }
  
  # -----------------------------------------------------------------------------
  # SORT
  ord <- order(meta$label)
  
  data <- data[ord, ]
  meta <- meta[ord, ]
  
  # -----------------------------------------------------------------------------
  # INDEX
  data.index = data[, index_names]
  data = data[, !(names(data) %in% index_names)]
  
  meta.index = meta[, index_names]
  meta = meta[, !(names(meta) %in% index_names)]
  
  stopifnot(all(meta.index == data.index))
  index = do.call(paste0, meta.index[index_names])
  
  # order = match(df1, rownames(df2))
  
  mat <- as.matrix(data)
  rownames(mat) <- index
  
  # -----------------------------------------------------------------------------
  # FILTER
  remove = meta$label == 'mix-vessels-PMN-MDSCs'
  mat = mat[!remove, ]
  meta = meta[!remove, ]
  
  # -----------------------------------------------------------------------------
  # AGGREGATION (mean per group)
  if (!is.null(aggregate_by)) {
    group <- meta[[aggregate_by]]
    mat <- rowsum(mat, group) / as.vector(table(group))
    meta <- meta[!duplicated(group), , drop = FALSE]
    stopifnot(all(meta[[aggregate_by]] == rownames(mat)))
  }
  
  # Define colour mappings for annotations
  label_colors <- unlist(colormaps$label)
  main_group_colors <- unlist(colormaps$main_group)
  patient_colors <- unlist(colormaps$pat_id)
  
  # -----------------------------------------------------------------------------
  # Create row annotations
  col_anno <- HeatmapAnnotation(
    label      = meta$label,
    main_group = meta$main_group,
    col = list(
      label      = label_colors,
      main_group = main_group_colors
    ),
    na_col = "#F0F0F0",
    show_legend = c(
      label = TRUE,
      main_group = TRUE
    )
  )
  
  # -----------------------------------------------------------------------------
  # Draw the heatmap
  if(column_split){
  heatmap_obj <- Heatmap(
    t(mat),   # transpose = rotate 90 degrees
    name = "Protein intensity",
    col = inferno(256),
    cluster_rows = TRUE,   # swap
    cluster_columns = cluster_rows,   # swap
    show_row_names = TRUE,            # was columns
    show_column_names = FALSE,        # was rows
    column_split = meta$main_group,   # was row_split
    column_title = NULL,              # was row_title
    top_annotation = col_anno,        # was right_annotation
    heatmap_legend_param = list(title = "Intensity")
  )}else{
    heatmap_obj <- Heatmap(
      t(mat),   # transpose = rotate 90 degrees
      name = "Protein intensity",
      col = inferno(256),
      cluster_rows = TRUE,   # swap
      cluster_columns = cluster_rows,   # swap
      show_row_names = TRUE,            # was columns
      show_column_names = FALSE,        # was rows
      # column_split = meta$main_group,   # was row_split
      column_title = NULL,              # was row_title
      top_annotation = col_anno,        # was right_annotation
      heatmap_legend_param = list(title = "Intensity"))
  }
  
  draw(heatmap_obj,
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
}

group_heatmap = function(main_group,
                         data_path,
                         metadata_path,
                         colormap_path,
                         index_names = c('sample_id', 'object_id'),
                         cluster_rows = FALSE,
                         subset_size = 2000,
                         replace_prefix = TRUE) {
  # -----------------------------------------------------------------------------
  # LOAD DATA
  data <- read_parquet(data_path)
  data = data |> select(-all_of(exclude_channels)) 
  
  meta <- read_parquet(metadata_path)
  colormaps <- yaml::read_yaml(colormap_path)
  
  # -----------------------------------------------------------------------------
  # SUBSET
  filter_ = meta$main_group == main_group
  meta = meta[filter_, ]
  data = data[filter_, ]
  
  if (nrow(data) > subset_size) {
    set.seed(0)
    subset = sample(nrow(data), subset_size)
    data = data[subset, ]
    meta = meta[subset, ]
  }
  
  # -----------------------------------------------------------------------------
  # SORT
  ord <- order(meta$label)
  
  data <- data[ord, ]
  meta <- meta[ord, ]
  
  # -----------------------------------------------------------------------------
  # INDEX
  data.index = data[, index_names]
  data = data[, !(names(data) %in% index_names)]
  
  meta.index = meta[, index_names]
  meta = meta[, !(names(meta) %in% index_names)]
  
  stopifnot(all(meta.index == data.index))
  index = do.call(paste0, meta.index[index_names])
  
  # order = match(df1, rownames(df2))
  
  mat <- as.matrix(data)
  rownames(mat) <- index
  
  # Define colour mappings for annotations
  label_colors <- unlist(colormaps$label)
  main_group_colors <- unlist(colormaps$main_group)
  patient_colors <- unlist(colormaps$pat_id)
  
  if(replace_prefix){
    pattern = paste0(main_group, "-")
    meta$label = str_replace_all(meta$label, pattern, "")    
    names(label_colors) <- str_replace_all(names(label_colors), pattern, "")
  }

  # -----------------------------------------------------------------------------
  # Create col annotations
  col_anno <- HeatmapAnnotation(
    label = meta$label,
    # main_group = meta$main_group,
    patient = meta$pat_id,
    col = list(
      label      = label_colors,
      # main_group = main_group_colors,
      patient    = patient_colors
    ),
    na_col = "#F0F0F0",
    show_legend = c(
      label = TRUE,
      # main_group = TRUE,
      patient = FALSE
    )
  )
  
  # -----------------------------------------------------------------------------
  # Draw the heatmap
  heatmap_obj <- Heatmap(
    t(mat),   # transpose = rotate 90 degrees
    name = "Protein intensity",
    col = inferno(256),
    cluster_rows = TRUE,   # swap
    cluster_columns = cluster_rows,   # swap
    show_row_names = TRUE,            # was columns
    show_column_names = FALSE,        # was rows
    column_split = meta$main_group,   # was row_split
    column_title = NULL,              # was row_title
    top_annotation = col_anno,        # was right_annotation
    heatmap_legend_param = list(title = "Intensity")
  )
  
  draw(heatmap_obj,
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
}


heatmap.caf = function(data_path,
                       metadata_path,
                       colormap_path,
                       index_names = c('sample_id', 'object_id'),
                       cluster_rows = TRUE) {
  # -----------------------------------------------------------------------------
  # LOAD DATA
  data <- read_parquet(data_path)
  data = data |> select(-all_of(exclude_channels)) 
  
  meta <- read_parquet(metadata_path)
  colormaps <- yaml::read_yaml(colormap_path)
  
  # -----------------------------------------------------------------------------
  # SORT
  ord <- order(meta$label)
  
  data <- data[ord, ]
  meta <- meta[ord, ]
  
  # -----------------------------------------------------------------------------
  # INDEX
  data.index = data[, index_names]
  data = data[, !(names(data) %in% index_names)]
  
  meta.index = meta[, index_names]
  meta = meta[, !(names(meta) %in% index_names)]
  
  stopifnot(all(meta.index == data.index))
  index = do.call(paste0, meta.index[index_names])
  
  # order = match(df1, rownames(df2))
  
  mat <- as.matrix(data)
  rownames(mat) <- index
  
  # -----------------------------------------------------------------------------
  # FILTER
  remove = meta$label == 'mix-vessels-PMN-MDSCs'
  mat = mat[!remove, ]
  meta = meta[!remove, ]
  
  filter1 = grepl("CAF", meta$label)
  # filter2 = grepl('pericytes', meta$label)
  # filter3 = grepl('stromal-', meta$label)
  filter_ = filter1
  mat = mat[filter_, ]
  meta = meta[filter_, ]
  
  caf_markers = c("smooth_muscle_actin",
                  "vimentin",
                  "collagen1",
                  "cd146",
                  "cnn1",
                  "cd105",
                  "ar",
                  "egr1",
                  "ces1"
                  )
  filter1 = colnames(mat) %in% caf_markers
  mat = mat[,filter1]
  mat = mat[,caf_markers]
  
  # -----------------------------------------------------------------------------
  # AGGREGATION (mean per group)
  aggregate_by = 'label'
  group <- meta[[aggregate_by]]
  mat <- rowsum(mat, group) / as.vector(table(group))
  meta <- meta[!duplicated(group), , drop = FALSE]
  stopifnot(all(meta[[aggregate_by]] == rownames(mat)))
  
  # Define colour mappings for annotations
  label_colors <- unlist(colormaps$label)
  main_group_colors <- unlist(colormaps$main_group)
  patient_colors <- unlist(colormaps$pat_id)
  
  pattern = paste0('stromal', "-")
  meta$label = str_replace_all(meta$label, pattern, "")    
  names(label_colors) <- str_replace_all(names(label_colors), pattern, "")
  
  # -----------------------------------------------------------------------------
  # Create row annotations
  col_anno <- HeatmapAnnotation(
    label      = meta$label,
    col = list(
      label = label_colors
    ),
    na_col = "#F0F0F0",
    show_legend = c(
      label = TRUE
    )
  )
  
  # -----------------------------------------------------------------------------
  # Draw the heatmap
  # mat = scale(mat)
  mat = scale(mat)
  heatmap_obj <- Heatmap(
    t(mat),   # transpose = rotate 90 degrees
    name = "Protein intensity",
    #col = inferno(256),
    col = circlize::colorRamp2(c(-2, 0, 2), c("lightblue", "white", "lightcoral")),
    cluster_rows = FALSE,   # swap
    cluster_columns = TRUE,   # swap
    show_row_names = TRUE,            # was columns
    show_column_names = FALSE,        # was rows
    # column_split = meta$main_group,   # was row_split
    column_title = NULL,              # was row_title
    top_annotation = col_anno,        # was right_annotation
    heatmap_legend_param = list(title = "Intensity"))
  
  draw(heatmap_obj,
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
}

# -----------------------------------------------------------------------------
# PLOT
ftype = '.pdf'

fname = paste0("heatmap-caf", ftype)
save_path = file.path(save_dir, fname)
# png(save_path, width = 10, height = 10, units = "in", res = 300)
pdf(save_path, width = 15, height = 10)
heatmap.caf(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  index_names = index_names,
  cluster_rows = FALSE
)
dev.off() 

fname = paste0("heatmap-cluster=false", ftype)
# png(save_path, width = 15, height = 10, units = "in", res = 300)
pdf(save_path, width = 15, height = 10)
save_path = file.path(save_dir, fname)
heatmap(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  index_names = index_names,
  cluster_rows = FALSE,
  subset = 2000
)
dev.off() 


fname = paste0("heatmap-cluster=true", ftype)
save_path = file.path(save_dir, fname)
# png(save_path, width = 15, height = 10, units = "in", res = 300)
pdf(save_path, width = 15, height = 10)
heatmap(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  index_names = index_names,
  cluster_rows = TRUE,
  subset = 2000
)
dev.off()  

fname = paste0("heatmap-cluster=true-agg=true", ftype)
save_path = file.path(save_dir, fname)
# png(save_path, width = 15, height = 10, units = "in", res = 300)
pdf(save_path, width = 15, height = 10)
heatmap.agg(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  index_names = index_names,
  cluster_rows = TRUE,
  aggregate_by = 'label',
  column_split = TRUE
)
dev.off() 

fname = paste0("heatmap-cluster=true-agg=true-split-false", ftype)
save_path = file.path(save_dir, fname)
# png(save_path, width = 15, height = 10, units = "in", res = 300)
pdf(save_path, width = 15, height = 10)
heatmap.agg(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  index_names = index_names,
  cluster_rows = TRUE,
  aggregate_by = 'label',
  column_split = FALSE
)
dev.off()  

groups = c('immune', 'epithelial', 'stromal', 'endothelial', 'undefined')
for(i in groups){
  fname = paste0("heatmap_", i, ".pdf")
  save_path = file.path(save_dir, fname)

  # png(save_path, width = 15, height = 10, units = "in", res = 300)
  pdf(save_path, width = 15, height = 10)
  group_heatmap(
    main_group = i,
    data_path = data_path,
    metadata_path = metadata_path,
    colormap_path = colormap_path,
    index_names = index_names,
    cluster_rows = FALSE,
    subset = 2000,
    replace_prefix = i != 'undefined'
  )
  dev.off()
}


groups = c('immune', 'epithelial', 'stromal', 'endothelial', 'undefined')
for(i in groups){
  fname = paste0("heatmap_", i, ftype)
  save_path = file.path(save_dir, fname)
  
  # png(save_path, width = 15, height = 10, units = "in", res = 300)
  pdf(save_path, width = 15, height = 10)
  group_heatmap(
    main_group = i,
    data_path = data_path,
    metadata_path = metadata_path,
    colormap_path = colormap_path,
    index_names = index_names,
    cluster_rows = FALSE,
    subset = 2000,
    replace_prefix = i != 'undefined'
  )
  dev.off()
}

# %% CAFs
ftype = '.pdf'
fname = paste0("heatmap-cafs-subset", ftype)
save_path = file.path(save_dir, fname)
pdf(save_path, width = 10, height = 10)
heatmap.caf(
  data_path = data_path,
  metadata_path = metadata_path,
  colormap_path = colormap_path,
  cluster_rows = TRUE,
)
dev.off()

