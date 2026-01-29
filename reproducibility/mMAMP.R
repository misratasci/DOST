
library(DOST)
set.seed(1999)

# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download mHypothalamus data from https://zenodo.org/records/10698931
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/"

# Helper function to load the data
load_mMAMP_sample <- function(dir.input, section = c("MA, MP")) {
  if (!(section %in% c("MA", "MP"))) {
    stop("section must be 'MA' or 'MP'")
  }
  dir.input <- file.path(dir.input, section)
  filename <- paste0(section, "_filtered_feature_bc_matrix.h5")
  sp_data <- Seurat::Load10X_Spatial(dir.input, filename = filename, filter.matrix = FALSE)
  df_meta <- read.table(file.path(dir.input, "metadata.tsv"),
                        sep = "\t", row.names = 1, header = TRUE)
  for (col in colnames(df_meta)) {
    sp_data <- SeuratObject::AddMetaData(sp_data,
                                         metadata = df_meta[col],
                                         col.name = col)
  }
  return(sp_data)
}

# Load data
section <- "MA"
sample <- load_mMAMP_sample(dir.input, section)
gt <- sample@meta.data$ground_truth
R <- length(unique(sample@meta.data$ground_truth))
X <- Seurat::GetAssayData(sample, layer = "counts")
coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]
results <- DOST(X, coords, R, refinement = FALSE)

# Plot results
plot_DOST_spatial(coords, results$labels)
plot_DOST_umap(results$Z, gt)

# Plot with Seurat
sample@meta.data$DOST <- results$labels
Seurat::SpatialDimPlot(sample, group.by = "DOST", pt.size.factor = 2.5, alpha = 1)

# Calculate ARI
mclust::adjustedRandIndex(results$labels, gt)
