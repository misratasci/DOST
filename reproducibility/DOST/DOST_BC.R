# DOST on the BC (human breast cancer, Visium) dataset.
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
#   BC_DOST_labels.RData       labels
#   BC_DOST_embeddings.RData   emb, the DOST embedding

library(DOST)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/BC"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"
# ---------------------------------------------------------------------------

set.seed(1999)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_BC_sample(dir.input)
gt <- sample@meta.data$fine_annot_type
R <- length(unique(gt))

X <- Seurat::GetAssayData(sample, layer = "counts")
coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]

# The BC domains are punctate rather than laminar, so no label refinement
results <- DOST(X, coords, R, refinement = FALSE)

ari <- mclust::adjustedRandIndex(results$labels, gt)
cat("Dataset: BC section1  ARI:", ari, "\n")

labels <- results$labels
save(labels, ari, file = file.path(dir.output, "BC_DOST_labels.RData"))
emb <- results$Z
save(emb, file = file.path(dir.output, "BC_DOST_embeddings.RData"))
