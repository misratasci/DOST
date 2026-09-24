# Runtime measurement for DOST on DLPFC12.
#
# Times the single DOST call, which covers normalization, gene selection,
# the distance matrices, the optimization and the clustering.
#
# Outputs written to dir.output:
#   DOST_runtime_slice<slice_index>.txt   one elapsed time in seconds per line,
#                                         appended (indexed 1-12, not by slice ID)

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Slice indices to time (1 to 12) and how many repeats per slice
slice_indices <- 1:12
n_repeats <- 5
# ---------------------------------------------------------------------------

set.seed(1999)

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

  for (i in 1:n_repeats) {
    t <- system.time(
      results <- DOST(X, coords, R)
    )
    cat("slice", slices[slice_index], "repeat", i, ":", t[3], "s\n")
    write(t[3],
          file = file.path(dir.output,
                           paste0("DOST_runtime_slice", slice_index, ".txt")),
          append = TRUE)
  }
}
