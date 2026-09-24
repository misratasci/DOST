# DOST on the DLPFC12 dataset.
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download DLPFC12 data from https://zenodo.org/records/10698880
#
# Expected layout:
#   DLPFC12/151673/{spatial/, 151673_filtered_feature_bc_matrix.h5,
#                   gt/tissue_positions_list_GTs.txt,
#                   gt/layered/151673_L1_barcodes.txt, ...}
#
# The DLPFC figure shows DOST with and without label refinement, so run this
# script twice, once with refinement TRUE and once with FALSE. The file name
# picks up a "_noref" tag in the second case.
#
# Outputs written to dir.output:
#   DOST_DLPFC.RData         results, aris   (refinement = TRUE)
#   DOST_DLPFC_noref.RData   results, aris   (refinement = FALSE)
#
# results is a list indexed by slice position 1-12, each entry the full DOST
# return value (Z, losses, labels). aris is one ARI per slice, in the same order.

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Slice indices to run (1 to 12).
slice_indices <- 1:12

# Label refinement
refinement <- TRUE
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

tag <- if (refinement) "" else "_noref"

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

results <- list()
aris <- c()

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

  results[[slice_index]] <- DOST(X, coords, R,
                                 refinement = refinement)

  ari <- mclust::adjustedRandIndex(results[[slice_index]]$labels,
                                   sample@meta.data$layers)
  aris <- c(aris, ari)

  cat("Slice:", slices[slice_index], " ARI:", ari, "\n")
}

save(results, aris,
     file = file.path(dir.output, paste0("DOST_DLPFC", tag, ".RData")))
