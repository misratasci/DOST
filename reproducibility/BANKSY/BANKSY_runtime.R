# Runtime measurement for BANKSY on DLPFC12.
#
# Times the BANKSY pipeline per slice: computeBanksy, runBanksyPCA, runBanksyUMAP and
# clusterBanksy. Loading, normalization and feature selection sit outside the timer.
#
# The timed clusterBanksy call uses the single resolution 0.7 from the Banksy vignette rather
# than the resolution sweep other scripts do.
#
# Outputs written to dir.output:
#   BANKSY_runtime<slice_index>.txt   one elapsed time in seconds per line, appended

library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(Seurat)
library(hdf5r)
source("reproducibility/BANKSY/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Slice indices to time (1 to 12) and how many repeats per slice
slice_indices <- 1:12
n_repeats <- 5
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

lambda <- c(0.2)
k_geom <- c(18, 18)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

for (slice_index in slice_indices) {

  sample <- load_DLPFC_sample(slices[slice_index], dir.input)

  gcm <- LayerData(sample, assay = "Spatial", layer = "counts")

  locs <- data.frame(sdimx = sample@meta.data$row, sdimy = sample@meta.data$col)
  row.names(locs) <- row.names(sample@meta.data)
  spatial_coor <- as.matrix(locs)
  rownames(spatial_coor) <- row.names(locs)

  se <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = spatial_coor)

  colData(se) <- DataFrame(
    sample_id = slices[slice_index],
    clust_annotation = factor(
      addNA(sample@meta.data$layers),
      exclude = NULL, labels = seq(length(unique(sample@meta.data$layers)))
    ),
    row.names = row.names(locs)
  )

  seu <- as.Seurat(se, data = NULL)
  seu <- NormalizeData(seu, scale.factor = 5000, normalization.method = "RC")
  feat <- VariableFeatures(FindVariableFeatures(seu, nfeatures = 2000))
  aname <- "normcounts"
  assay(se, aname) <- GetAssayData(seu)
  se <- se[feat, ]
  rm(seu)

  se0 <- se

  for (i in 1:n_repeats) {
    se <- se0
    set.seed(1000)
    t <- system.time({
      se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
      se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
      se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
      se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 0.7)
    })
    cat("slice", slices[slice_index], "repeat", i, ":", t[3], "s\n")
    write(t[3],
          file = file.path(dir.output, paste0("BANKSY_runtime", slice_index, ".txt")),
          append = TRUE)
  }
}
