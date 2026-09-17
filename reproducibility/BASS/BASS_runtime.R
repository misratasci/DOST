# Runtime measurement for BASS on DLPFC12.
#
# Times the full pipeline per slice:
# createBASSObject, BASS.preprocess, BASS.run and BASS.postprocess.
# Data loading sits outside the timer.
#
# Outputs written to dir.output:
#   BASS_runtime<slice_index>.txt   one elapsed time in seconds per line, appended
#                                   (indexed 1-12, not by slice ID)

library(BASS)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"

# Slice indices to time (1 to 12) and how many repeats per slice
slice_indices <- 1:12
n_repeats <- 5
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

for (slice_index in slice_indices) {

  slice_nos <- c(slice_index)
  samples <- lapply(slice_nos, function(i) {
    load_DLPFC_sample(slices[i], dir.input)
  })
  cnts <- lapply(samples, function(i) {
    i@assays$Spatial$counts
  })
  xym <- lapply(samples, function(i) {
    data.frame(row = i@meta.data$row, col = i@meta.data$col)
  })
  for (i in 1:length(slice_nos)) {
    row.names(xym[[i]]) <- row.names(samples[[i]]@meta.data)
  }

  C <- 20
  if (slice_index > 4 & slice_index < 9) {
    R <- 5
  } else {
    R <- 7
  }

  for (i in 1:n_repeats) {
    set.seed(1999)
    t <- system.time({
      BASS <- createBASSObject(cnts, xym, C = C, R = R)
      BASS <- BASS.preprocess(BASS)
      BASS <- BASS.run(BASS)
      BASS <- BASS.postprocess(BASS)
    })
    cat("slice", slices[slice_index], "repeat", i, ":", t[3], "s\n")
    write(t[3],
          file = file.path(dir.output, paste0("BASS_runtime", slice_index, ".txt")),
          append = TRUE)
  }
}
