# Ablation: the spatial weight lambda across its whole range, on DLPFC12.
#
# Two passes over all 12 slices. The first sweeps lambda with refinement off and
# saves the raw labels. The second loads those, applies label refinement to each
# lambda's labels and saves them again
#
# Outputs written to dir.output:
#   DOST_DLPFC<slice_index>_lambda_noref.RData   results, aris
#   DOST_DLPFC<slice_index>_lambda_ref.RData     results, aris
#
# results is a list keyed by lambda as a string, each entry the full DOST return
# value (Z, losses, labels). aris is named the same way. Note both are indexed
# by slice position 1-12, not by slice ID.

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Lambda grid to sweep
lambdas <- seq(0, 1, by = 0.1)

# Slice indices to run (1 to 12)
slice_indices <- 1:12
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

domains_for <- function(slice_index) {
  if (slice_index > 4 & slice_index < 9) 5 else 7
}

# ---------------------------------------------------------------------------
# Pass 1: sweep lambda without refinement
# ---------------------------------------------------------------------------

for (slice_index in slice_indices) {

  cat(paste0("\nSlice ", slice_index, "\n"))

  R <- domains_for(slice_index)

  sample <- load_DLPFC_sample(slices[slice_index], dir.input)
  X <- Seurat::GetAssayData(sample, layer = "counts")
  coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]

  results <- list()
  aris <- c()

  for (lambda in lambdas) {
    results[[as.character(lambda)]] <-
      DOST(X, coords, R,
           lambda = lambda,
           refinement = FALSE)

    aris[[as.character(lambda)]] <- mclust::adjustedRandIndex(
      results[[as.character(lambda)]]$labels, sample@meta.data$layers)

    cat("slice:", slices[slice_index], " lambda:", lambda,
        " ARI:", aris[[as.character(lambda)]], "\n")
  }

  save(results, aris,
       file = file.path(dir.output,
                        paste0("DOST_DLPFC", slice_index, "_lambda_noref.RData")))
}

# ---------------------------------------------------------------------------
# Pass 2: refine the labels from pass 1
# ---------------------------------------------------------------------------

for (slice_index in slice_indices) {

  cat(paste0("\nSlice ", slice_index, " (refining)\n"))

  sample <- load_DLPFC_sample(slices[slice_index], dir.input)
  coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]

  env <- new.env()
  load(file.path(dir.output,
                 paste0("DOST_DLPFC", slice_index, "_lambda_noref.RData")),
       envir = env)
  results <- env$results
  aris <- env$aris

  for (lambda in lambdas) {
    results[[as.character(lambda)]]$labels <-
      DOST:::refine_labels(coords, results[[as.character(lambda)]]$labels)

    aris[[as.character(lambda)]] <- mclust::adjustedRandIndex(
      results[[as.character(lambda)]]$labels, sample@meta.data$layers)

    cat("slice:", slices[slice_index], " lambda:", lambda,
        " refined ARI:", aris[[as.character(lambda)]], "\n")
  }

  save(results, aris,
       file = file.path(dir.output,
                        paste0("DOST_DLPFC", slice_index, "_lambda_ref.RData")))
}
