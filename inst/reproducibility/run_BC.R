
library(DOST)
set.seed(1999)

# Change this path to where you downloaded the data
dir.input <- "path/to/data/"

# Helper function to load the data
load_BC_sample <- function(dir.input, section_no = 1) {
  if (section_no != 1 & section_no != 2) {
    stop("section_no must be 1 or 2")
  }
  filename <- paste0("section", section_no, "_filtered_feature_bc_matrix.h5")
  sp_data <- Seurat::Load10X_Spatial(file.path(dir.input, paste0("section", section_no)), filename = filename, filter.matrix = FALSE)
  if (section_no == 1) {
    df_meta <- read.table(file.path(dir.input, paste0("section", section_no), "gt", "gold_metadata.tsv"),
                          sep = "\t", row.names = 1, header = TRUE)
    sp_data <- SeuratObject::AddMetaData(sp_data,
                                         metadata = df_meta[1],
                                         col.name = 'annot_type')
    sp_data <- SeuratObject::AddMetaData(sp_data,
                                         metadata = df_meta[2],
                                         col.name = 'fine_annot_type')
  }
  return(sp_data)
}

gt <- sample@meta.data$fine_annot_type
sample <- load_BC_sample(dir.input)
R <- length(unique(gt))
X <- Seurat::GetAssayData(sample, layer = "counts")
coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]
results <- DOST(X, coords, R, lambda = 0.03, lr = 16, refinement = FALSE)

plot_DOST_spatial(coords, results$labels)
plot_DOST_umap(results$Z, gt)

mclust::adjustedRandIndex(results$labels, gt)
