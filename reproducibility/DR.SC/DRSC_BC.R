# DR.SC on the BC dataset.
#
# Method: https://feiyoung.github.io/DR.SC/
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download BC data from https://zenodo.org/records/10698903
#
# Expected layout:
#   BC/section1/{spatial/, section1_filtered_feature_bc_matrix.h5,
#                gt/gold_metadata.tsv, gt/tissue_positions_list_GTs.txt}
#
# Outputs written to dir.output:
#   BC_DRSC_labels.RData       labels, the spatial.drsc.cluster assignment
#   BC_DRSC_embeddings.RData   emb, the "dr-sc" reduction

library(Seurat)
library(DR.SC)
source("reproducibility/DR.SC/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/BC"

# Change this path to where you want to save the results
dir.output <- "path/to/output"
# ---------------------------------------------------------------------------

# Number of highly variable genes DR.SC is given
n_features <- 500

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_BC_sample(dir.input)

cluster.number <- length(unique(sample@meta.data$fine_annot_type))

# DR.SC expects the spatial coordinates to be in the Seurat object metadata, so we add them here
df_meta <- read.table(file.path(dir.input, "section1", "gt", "tissue_positions_list_GTs.txt"),
                      sep = ",", row.names = 1)
xym <- data.frame(row = df_meta['V3'], col = df_meta['V4'])
sample@meta.data$row <- xym$V3
sample@meta.data$col <- xym$V4

sample <- NormalizeData(sample, verbose = FALSE)
seu <- FindVariableFeatures(sample, nfeatures = n_features, verbose = FALSE)

seu <- DR.SC(seu, K = as.numeric(cluster.number), platform = 'Visium', verbose = FALSE)

ari_drsc <- mclust::adjustedRandIndex(seu$spatial.drsc.cluster, seu@meta.data$fine_annot_type)
print(paste0("ARI: ", ari_drsc))

emb <- Embeddings(seu, reduction = "dr-sc")
save(emb, file = file.path(dir.output,
                           paste0("BC_DRSC_embeddings.RData")))
labels <- seu$spatial.drsc.cluster
save(labels, file = file.path(dir.output,
                              paste0("BC_DRSC_labels.RData")))
