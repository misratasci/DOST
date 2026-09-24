# Ablation: how long DOST needs to optimize, on one DLPFC12 slice.
#
# Outputs written to dir.output:
#   DOST_DLPFC_ablation_max_iter_<slice_index>.RData   results
#
# results is a list keyed by the iteration budget as a string, each entry the
# full DOST return value (Z, losses, labels).

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Iteration budgets to sweep
max_iter_grid <- seq(0, 60, by = 20)

# Slice to run (1 to 12)
slice_index <- 9
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

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

for (max_iterations in max_iter_grid) {
  results[[as.character(max_iterations)]] <-
    DOST(X, coords, R,
         max_iterations = max_iterations)

  ari <- mclust::adjustedRandIndex(results[[as.character(max_iterations)]]$labels,
                                   sample@meta.data$layers)
  cat("slice:", slices[slice_index], " max_iterations:", max_iterations,
      " ARI:", ari, "\n")
}

save(results,
     file = file.path(dir.output,
                      paste0("DOST_DLPFC_ablation_max_iter_", slice_index,".RData")))
