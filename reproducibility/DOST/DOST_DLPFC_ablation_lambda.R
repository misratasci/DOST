# Ablation: the spatial weight lambda around its default, on DLPFC12.
#
# Outputs written to dir.output:
#   DOST_DLPFC_ablation_lambda_<slice_index>.RData   results
#
# results is a list keyed by lambda as a string, each entry the full DOST return
# value (Z, losses, labels).

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Lambda grid to sweep
lambdas <- seq(0, 0.09, by = 0.03)

# Slices to run (1 to 12), one output file each
# Slice 1 and 9 were run for the ablation study
slice_indices <- c(1, 9)
# ---------------------------------------------------------------------------

# Slice IDs
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

  for (lambda in lambdas) {
    results[[as.character(lambda)]] <-
      DOST(X, coords, R,
           lambda = lambda)

    ari <- mclust::adjustedRandIndex(results[[as.character(lambda)]]$labels,
                                     sample@meta.data$layers)
    cat("slice:", slices[slice_index], " lambda:", lambda, " ARI:", ari, "\n")
  }

  save(results,
       file = file.path(dir.output,
                        paste0("DOST_DLPFC_ablation_lambda_", slice_index,
                               ".RData")))
}
