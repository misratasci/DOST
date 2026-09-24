# Ablation: the size of the DOST embedding, on DLPFC12.
#
# Outputs written to dir.output:
#   DOST_DLPFC<slice_index>_embdim.RData   results, aris
#
# results is a list keyed by the embedding dimension as a string, each entry
# the full DOST return value (Z, losses, labels). aris is named the same way.
# Note both are indexed by slice position 1-12, not by slice ID.

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Embedding dimensions to sweep
embdims <- c(10, 20, 30, 40)

# Slice indices to run (1 to 12)
slice_indices <- 1:12
# ---------------------------------------------------------------------------

slices <- c(151507:151510, 151669:151676)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

for (slice_index in slice_indices) {

  cat(paste0("\nSlice ", slice_index, "\n"))

  if (slice_index > 4 & slice_index < 9) {
    R <- 5
  } else {
    R <- 7
  }

  sample <- load_DLPFC_sample(slices[slice_index], dir.input)
  X <- Seurat::GetAssayData(sample, layer = "counts")
  coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]

  results <- list()
  aris <- c()

  for (embedding_dim in embdims) {
    results[[as.character(embedding_dim)]] <-
      DOST(X, coords, R,
           embedding_dim = embedding_dim)

    aris[[as.character(embedding_dim)]] <- mclust::adjustedRandIndex(
      results[[as.character(embedding_dim)]]$labels, sample@meta.data$layers)

    cat("slice:", slices[slice_index], " embedding_dim:", embedding_dim,
        " ARI:", aris[[as.character(embedding_dim)]], "\n")
  }

  save(results, aris,
       file = file.path(dir.output,
                        paste0("DOST_DLPFC", slice_index, "_embdim.RData")))
}
