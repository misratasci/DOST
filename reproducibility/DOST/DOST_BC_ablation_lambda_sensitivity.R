# Ablation: the spatial weight lambda across its whole range, on BC.
#
# Sweeps lambda over the grid below with refinement off, matching the BC
# setting in DOST_BC.R, and scores each run against the fine annotation.
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
#   DOST_BC_lambda_sensitivity.RData   results, aris
#
# results is a list keyed by lambda as a string, each entry the full DOST return
# value (Z, losses, labels). aris is named the same way.

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/BC"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Lambda grid to sweep
lambdas <- seq(0, 0.1, by = 0.01)
# ---------------------------------------------------------------------------

set.seed(1999)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_BC_sample(dir.input)
gt <- sample@meta.data$fine_annot_type
R <- length(unique(gt))

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
    results[[as.character(lambda)]]$labels, gt)

  cat("Dataset: BC section1  lambda:", lambda,
      " ARI:", aris[[as.character(lambda)]], "\n")
}

save(results, aris,
     file = file.path(dir.output, "DOST_BC_lambda_sensitivity.RData"))
