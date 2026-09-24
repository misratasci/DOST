# DOST on the mMAMP (mouse brain anterior, "MA") dataset.
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
#   mMAMP_DOST_labels.RData       labels
#   mMAMP_DOST_embeddings.RData   emb, the DOST embedding

library(DOST)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mMAMP"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"
# ---------------------------------------------------------------------------

set.seed(1999)

section <- "MA"

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_mMAMP_sample(dir.input, section)
gt <- sample@meta.data$ground_truth
R <- length(unique(gt))

X <- Seurat::GetAssayData(sample, layer = "counts")
coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]

# The mMAMP domains are punctate rather than laminar, so no label refinement
results <- DOST(X, coords, R, refinement = FALSE)

ari <- mclust::adjustedRandIndex(results$labels, gt)
cat("Dataset: mMAMP", section, " ARI:", ari, "\n")

labels <- results$labels
save(labels, ari, file = file.path(dir.output, "mMAMP_DOST_labels.RData"))
emb <- results$Z
save(emb, file = file.path(dir.output, "mMAMP_DOST_embeddings.RData"))
