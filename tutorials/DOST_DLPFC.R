
library(DOST)
set.seed(1999)

# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download DLPFC12 data from https://zenodo.org/records/10698880
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/"

# Helper function to load the data
load_DLPFC_sample <- function(slice.id, dir.input) {
  filename <- paste0(slice.id, "_filtered_feature_bc_matrix.h5")
  data.dir <- file.path(dir.input, slice.id)
  sp_data <- Seurat::Load10X_Spatial(data.dir, filename = filename, filter.matrix = FALSE)
  #add the annotations
  df_meta <- read.table(file.path(data.dir, 'gt', 'tissue_positions_list_GTs.txt'),
                        sep = ",", row.names = 1)
  common_cells <- colnames(sp_data[["Spatial"]]) %in% rownames(df_meta)
  sp_data <- sp_data[, common_cells]
  layer.data <- data.frame()
  layers <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM')
  for (l in layers) {
    filename <- paste0(slice.id, "_", l, "_barcodes.txt")
    filename <- file.path(data.dir, 'gt', 'layered', filename)
    if (!file.exists(filename)) next
    data.temp <- read.table(filename)
    data.temp <- data.frame(barcode = data.temp[,1], layer = l, row.names = data.temp[,1])
    layer.data <- rbind(layer.data, data.temp)
  }
  sp_data <- SeuratObject::AddMetaData(sp_data,
                                       metadata = df_meta['V3'],
                                       col.name = 'row')
  sp_data <- SeuratObject::AddMetaData(sp_data,
                                       metadata = df_meta['V4'],
                                       col.name = 'col')
  sp_data <- SeuratObject::AddMetaData(sp_data,
                                       metadata = layer.data['layer'],
                                       col.name = 'layers')
  return (sp_data)
}

# Slice IDs
slices <- c(151507:151510, 151669:151676)

# Run DOST for a selected slice index (from 1 to 12)
slice_index <- 9

# Load data
sample <- load_DLPFC_sample(slices[slice_index], dir.input)
gt <- sample@meta.data$layers
R <- length(unique(gt))
X <- Seurat::GetAssayData(sample, layer = "counts")
coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]

# Run DOST
results <- DOST(X, coords, R, refinement = TRUE)

# Plot results
plot_DOST_spatial(coords, results$labels)
plot_DOST_umap(results$Z, gt)

# Plot with Seurat
sample@meta.data$DOST <- results$labels
Seurat::SpatialDimPlot(sample, group.by = "DOST", pt.size.factor = 2.5, alpha = 1)

# Calculate ARI
mclust::adjustedRandIndex(results$labels, gt)
