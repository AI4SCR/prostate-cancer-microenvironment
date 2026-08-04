library("imcRtools")
library("CATALYST")
library("pheatmap")
library("dittoSeq")
library("patchwork")
library("cytomapper")
library("stringr")
library('assert')
library('tidyverse')

base_dir = '/home/labadmin/data/pca-v3'
# base_dir = '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/remotes/unibe/data/pca-v3'

spillover_acquisition = 'spillslide_pca_panel_190224'
path.spillover.dir = file.path(base_dir, 'spillover', spillover_acquisition)
path.spillover.matrix = file.path(path.spillover.dir, 'spillover_matrix.RDS')
path.panel = file.path(base_dir, 'images', 'raw', 'panel.csv')
prefix = file.path(base_dir, "report", "spillover", spillover_acquisition)

path.analytics = file.path(prefix)

if (!dir.exists(path.analytics)) {
  dir.create(path.analytics, showWarnings = FALSE, recursive = TRUE)
}

panel = read_csv(path.panel)

# NOTE: rows are the measured channels, columns are measured pixles. `colData(sce)$sample_id` indicates from which spot the pixel comes
sce <- readSCEfromTXT(path.spillover.dir, read_metal_from_filename = F)

pattern = "_([A-z]{2,3})(\\d{2,3})_\\d+$"
sce$sample_metal <- str_to_title(str_match(sce$sample_id, pattern)[,2])
assert(sum(is.na(sce$sample_metal)) == 0)

sce$sample_mass <- str_match(sce$sample_id, pattern)[,3]
assert(sum(is.na(sce$sample_mass)) == 0)

sce$sample_id <- paste0(sce$sample_metal, sce$sample_mass)
assert(sum(is.na(sce$sample_id)) == 0)

indices <- sub(".*\\.", "", rownames(colData(sce)))
rownames(colData(sce)) <- paste(sce$sample_id, indices, sep = ".")

pattern = "([A-z]{2,3})(\\d{2,3})"
rowData(sce)$metal_mass = str_match(rowData(sce)$marker_name, pattern)[,3]
rowData(sce)$metal_name = str_match(rowData(sce)$marker_name, pattern)[,2]

# NOTE: we only keep the metal masses that were stained
assert(all(panel$metal_mass %in% rowData(sce)$metal_mass) == T)
mask = rowData(sce)$metal_mass %in% panel[panel$is_stained,]$metal_mass
sce = sce[mask,]

# NOTE: we use cofactor 1 and not 5 like in the Bodenmiller tutorial
cofactor = 1
assay(sce, "exprs") <- asinh(counts(sce)/cofactor)


# NOTE: rows are prepared and measured isotopes, columns are measured channels
pdf(file=file.path(path.analytics, "spot-heatmap.pdf"))
plotSpotHeatmap(sce)
dev.off()

pdf(file=file.path(path.analytics, "spot-heatmap-thres-200.pdf"))
plotSpotHeatmap(sce, log = F, threshold = 200)
dev.off()

pdf(file=file.path(path.analytics, "spot-heatmap-thres-100.pdf"))
plotSpotHeatmap(sce, log = F, threshold = 100)
dev.off()

# De-barcoding
bc_key <- as.numeric(unique(sce$sample_mass))
bc_key <- bc_key[order(bc_key)]

sce <- assignPrelim(sce, bc_key = bc_key)
sce <- estCutoffs(sce)
sce <- applyCutoffs(sce)

cur_table <- table(sce$bc_id, sce$sample_mass)
cur_table[1:10, 1:10]

# Visualize the correctly and incorrectly assigned pixels
pdf(file=file.path(path.analytics, 'pixel-assignment.pdf'))
pheatmap(log10(cur_table + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

# Compute the fraction of unassigned pixels per spot
cur_frac = cur_table["0",] / colSums(cur_table)
cur_frac
print(c(min(cur_frac), max(cur_frac)))

# NOTE: if very few pixels (~100) are measured per spot, the minevents parameter value needs to be lowered.
# sce <- filterPixels(sce, minevents = 40, correct_pixels = T)

# NOTE: docs say: Returns a square compensation matrix with dimensions and dimension names matching those of the input flowFrame. But dim(sce) == 49x70000, dim(sm) = 35 x 49
sce <- computeSpillmat(sce)

pdf(file=file.path(path.analytics, 'spillmat.pdf'))
plotSpillmat(sce)
dev.off()

sm <- metadata(sce)$spillover_matrix
saveRDS(sm, file = path.spillover.matrix)

ggdata = sm %>% 
  as.data.frame %>% 
  rownames_to_column("source") %>% 
  pivot_longer(-c(source), names_to = 'target')
ggdata$value[ggdata$target == ggdata$source] = NA

pattern = "([A-za-z]{2,3})(\\d{2,3})"
ggdata$source_mass = str_match(ggdata$source, pattern)[,3]
ggdata$target_mass = str_match(ggdata$target, pattern)[,3]

panel = panel %>% arrange(metal_mass)
metal2name = setNames(panel$name, panel$metal_mass)
ggdata$source_name = factor(
  metal2name[as.character(ggdata$source_mass)],
  levels = rev(panel$name)
  )
ggdata$target_name = factor(
  metal2name[as.character(ggdata$target_mass)],
  levels = panel$name
)

saveRDS(ggdata, file=file.path(prefix, 'ggdata.RDS'))

ggplot(ggdata, aes(target_name, source_name, fill= value)) + 
  geom_tile(color = "black") + 
  geom_text(data = subset(ggdata, value > 0), aes(label = sprintf("%0.2f", value * 100)), size = 2) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(path.analytics, "spillover-custom.pdf"), device = "pdf", width = 10, height = 8)

pdf(file.path(path.analytics, 'hist.pdf'), width = 14, height = 7)
hist(ggdata$value[ggdata$value > 0])
dev.off()
