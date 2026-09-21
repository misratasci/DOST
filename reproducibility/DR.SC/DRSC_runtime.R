# Runtime measurement for DR.SC on DLPFC12.
#
# Counterpart to BASS_runtime.R, BANKSY_runtime.R and the Python methods' runtime
# scripts. Times normalization, variable-feature selection and the DR.SC fit. Data
# loading sits outside the timer, matching the original script.
#
# Outputs written to dir.output:
#   DRSC_runtime<slice_index>.txt   one elapsed time in seconds per line, appended
#                                   (indexed 1-12, not by slice ID)

library(Seurat)
library(DR.SC)
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

n_features <- 500

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

for (slice_index in slice_indices) {

  if (slice_index > 4 & slice_index < 9) {
    cluster.number <- 5
  } else {
    cluster.number <- 7
  }

  sample <- load_DLPFC_sample(slices[slice_index], dir.input)

  for (i in 1:n_repeats) {
    t <- system.time({
      sample_i <- NormalizeData(sample, verbose = FALSE)
      seu <- FindVariableFeatures(sample_i, nfeatures = n_features, verbose = FALSE)
      seu <- DR.SC(seu, K = as.numeric(cluster.number), platform = 'Visium', verbose = FALSE)
    })
    cat("slice", slices[slice_index], "repeat", i, ":", t[3], "s\n")
    write(t[3],
          file = file.path(dir.output, paste0("DRSC_runtime", slice_index, ".txt")),
          append = TRUE)
  }
}
