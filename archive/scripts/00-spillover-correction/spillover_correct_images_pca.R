# nohup Rscript ~/projects/BLCaPCa/create-report/spillover_correct_images_pca.R > ~/data/pca/logs/spillover_correct_images_pca.out 2>&1 &

library(CATALYST)
library(cytomapper)
library(BiocParallel)
library(stringr)
library(tidyverse)
library(ijtiff)
library(lubridate)

# Define the function
get_name_of_closer_acquisition <- function(filename, dates) {
  # Extract the date string from the filename
  date_str <- sub("_(.*)", "", filename)
  
  # Convert the extracted date string to a Date object
  file_date <- ymd(date_str)
  
  # Compute the absolute difference between the file date and each of the other dates
  diff <- abs(as.numeric(difftime(file_date, dates, units = "days")))
  return(names(dates)[which.min(diff)])
}


base_dir = '/home/labadmin/data/pca-v3'

spillover_acquisition = c('spillslide_pca_panel_120124','spillslide_pca_panel_190224')
spillover_acquisition.dates = spillover_acquisition %>% str_split('_') %>% dmy
names(spillover_acquisition.dates) = spillover_acquisition

path.spillover.dir = file.path(base_dir, 'spillover', spillover_acquisition)
path.spillover.matrix = file.path(path.spillover.dir, 'spillover_matrix.RDS')

sm = lapply(path.spillover.matrix, readRDS)
names(sm) = spillover_acquisition

path.img.dir = file.path(base_dir, 'images', 'filtered')

path.img.dir.compensated = file.path(base_dir, 'images', 'compensated')
if (!dir.exists(path.img.dir.compensated)) {
  dir.create(path.img.dir.compensated, showWarnings = FALSE, recursive = TRUE)
}

path.panel = file.path(base_dir, 'images', 'filtered', 'panel.csv')
panel = read_csv(path.panel)
panel$name_for_compensation = paste0(panel$metal_name, panel$metal_mass, 'Di')
channel.names.stained = panel[panel$is_stained,]$name_for_compensation

# note: the missing interactions are from DNA and ICSK channels
sm.adapted = lapply(sm, function(x) adaptSpillmat(x, channel.names.stained))

img.filenames = list.files(path.img.dir, pattern = '.tiff')
N = length(img.filenames)
i = 1
for(img.name in img.filenames){
  # img.name = img.filenames[1]
  print(paste('(', i, '/', N, ')', img.name))

  img = loadImages(path.img.dir, pattern = img.name)
  channelNames(img) = panel$name_for_compensation
  img.stained = getChannels(img, i = channel.names.stained)

  spillover_acquisition.name = get_name_of_closer_acquisition(img.name, spillover_acquisition.dates)
  sm.current = sm.adapted[[spillover_acquisition.name]]
  img.compensated <- compImage(img.stained, sm.current, 
                           BPPARAM = MulticoreParam())
  
  sample.name = str_replace(img.name, '.tiff', '')
  img.data = imageData(img.compensated[[sample.name]])
  img.path = file.path(path.img.dir.compensated, img.name)
  img.dim = dim(img.data)
  img.ijtiff = as_ijtiff_img(img.data, dim=c(img.dim,1)) # note: we add a forth plane dimension
  # img.ijtiff = as_ijtiff_img(array(as.integer(img.data), dim = img.dim), dim=c(img.dim,1))
  
  write_tif(img.ijtiff, img.path, overwrite=T)
  i = i + 1
}

file.copy(file.path(path.img.dir, 'panel.csv'), file.path(path.img.dir.compensated, 'panel.csv'))
file.copy(file.path(path.img.dir, 'images.csv'), file.path(path.img.dir.compensated, 'images.csv'))
