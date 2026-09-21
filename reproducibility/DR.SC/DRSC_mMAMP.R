# DR.SC on the mMAMP dataset.
#
# Method: https://feiyoung.github.io/DR.SC/
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download mMAMP data from https://zenodo.org/records/10698931
#
# Expected layout:
#   mMAMP/MA/{spatial/, MA_filtered_feature_bc_matrix.h5, metadata.tsv,
#             gt/tissue_positions_list_GTs.txt}
#
# Outputs written to dir.output:
#   mMAMP_DRSC_labels.RData       labels, the spatial.drsc.cluster assignment
#   mMAMP_DRSC_embeddings.RData   emb, the "dr-sc" reduction

library(Seurat)
library(DR.SC)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mMAMP"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"
# ---------------------------------------------------------------------------

# Number of highly variable genes DR.SC is given
n_features <- 500

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_mMAMP_sample(dir.input)

cluster.number <- length(unique(sample@meta.data$ground_truth))

# DR.SC expects the spatial coordinates to be in the Seurat object metadata, so we add them here
df_meta <- read.table(file.path(dir.input, "MA", "gt", "tissue_positions_list_GTs.txt"),
                      sep = "\t", header = TRUE)
sample@meta.data$row <- df_meta$array_row
sample@meta.data$col <- df_meta$array_col

sample <- NormalizeData(sample, verbose = FALSE)
seu <- FindVariableFeatures(sample, nfeatures = n_features, verbose = FALSE)

seu <- DR.SC(seu, K = as.numeric(cluster.number), platform = 'Visium', verbose = FALSE)

ari_drsc <- mclust::adjustedRandIndex(seu$spatial.drsc.cluster, seu@meta.data$ground_truth)
print(paste0("ARI: ", ari_drsc))

emb <- Embeddings(seu, reduction = "dr-sc")
save(emb, file = file.path(dir.output,
                           paste0("mMAMP_DRSC_embeddings.RData")))
labels <- seu$spatial.drsc.cluster
save(labels, file = file.path(dir.output,
                              paste0("mMAMP_DRSC_labels.RData")))
