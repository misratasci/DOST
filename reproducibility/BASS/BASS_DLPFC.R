# BASS on the DLPFC12 dataset.
#
# Adapted from https://benchmarkst-reproducibility.readthedocs.io/en/latest/BASS_clustering.html
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download DLPFC12 data from https://zenodo.org/records/10698880
#
# Outputs written to dir.output:
#   DLPFC_<slice>_BASS_labels.RData   zlabels (domains) and clabels (cell types)
#   bass_aris.RData                   bass_aris, one ARI per slice covered by the loop

library(BASS)
source("reproducibility/BASS/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Slice indices to run (1 to 12). Set to a single index, e.g. 9, for one slice;
# leave as 1:12 to produce the full bass_aris vector
slice_indices <- 1:12
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

bass_aris <- c()

for (slice_index in slice_indices) {

  set.seed(1999)

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

  # C is the number of cell types, R the number of spatial domains
  C <- 20
  if (slice_index > 4 & slice_index < 9) {
    R <- 5
  } else {
    R <- 7
  }

  BASS <- createBASSObject(cnts, xym, C = C, R = R,
                           beta_method = "SW",
                           init_method = "mclust",
                           nsample = 10000)
  listAllHyper(BASS)
  BASS <- BASS.preprocess(BASS)
  BASS <- BASS.run(BASS)
  BASS <- BASS.postprocess(BASS)

  # Save the labels
  zlabels <- BASS@results$z
  clabels <- BASS@results$c
  save(zlabels, clabels,
       file = file.path(dir.output,
                        paste0("DLPFC_", slices[slice_index], "_BASS_labels.RData")))

  gt <- samples[[1]]@meta.data$layers
  ari <- mclust::adjustedRandIndex(zlabels[[1]], gt)
  bass_aris <- c(bass_aris, ari)

  cat("Slice:", slices[slice_index], " ARI:", ari, "\n")
}

save(bass_aris, file = file.path(dir.output, "bass_aris.RData"))
