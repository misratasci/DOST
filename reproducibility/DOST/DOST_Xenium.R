# DOST on the Xenium Human Breast Gene Expression test dataset.
#
# Downloaded from the official 10x Genomics website:
# curl -O https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip
#
# Outputs written to dir.output:
# - Xenium_DOST_labels.RData: cluster labels for each cell
# - Xenium_DOST_embeddings.RData: low-dimensional embedding for each cell

library(Seurat)
library(DOST)

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/Xenium_V1_human_Breast_2fov_outs"

# Change this path to where you want to save the results
dir.output <- "path/to/output"
# ---------------------------------------------------------------------------

set.seed(1999)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- LoadXenium(data.dir = dir.input)

R <- 5 # Chosen as an exploratory value

X <- GetAssayData(sample, assay = "Xenium", layer = "counts")
coords <- GetTissueCoordinates(sample)

# Remove cells with zero counts
cell_sums <- colSums(X)
keep <- cell_sums > 0
sample <- subset(sample, cells = colnames(X)[keep])
X <- X[, keep]
coords <- coords[keep, ]
coords <- coords[, c("x", "y")]

results <- DOST(
  X = X,
  coords = coords,
  R = R,
  selected_genes = "all", # Only 280 genes are detected in the panel
  refinement = FALSE
)

labels <- results$labels
save(labels, file = file.path(dir.output, "Xenium_DOST_labels.RData"))
emb <- results$Z
save(emb, file = file.path(dir.output, "Xenium_DOST_embeddings.RData"))
