# BANKSY on the DLPFC12 dataset.
#
# Method: https://prabhakarlab.github.io/Banksy/
# Preprocessing follows https://prabhakarlab.github.io/Banksy/articles/multi-sample.html
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download DLPFC12 data from https://zenodo.org/records/10698880
#
# Outputs written to dir.output:
#   DLPFC_BANKSY_slice_<slice>_labels.RData       labels, the chosen clustering
#   DLPFC_BANKSY_slice_<slice>_embeddings.RData   emb, the BANKSY PCA embedding
#   DLPFC_BANKSY_aris.RData                       banksy_aris, one ARI per slice in the loop

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

# Slice indices to run (1 to 12). Leave as 1:12 to produce the full banksy_aris
# vector the DLPFC figure expects.
slice_indices <- 1:12
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

# BANKSY recommends the following hyperparameters for 10x Visium v1v2
# https://prabhakarlab.github.io/Banksy/articles/parameter-selection.html
lambda <- c(0.2)
k_geom <- c(18, 18)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

banksy_aris <- c()

for (slice_index in slice_indices) {

  if (slice_index > 4 & slice_index < 9) {
    cluster.number <- 5
  } else {
    cluster.number <- 7
  }

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

  se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

  set.seed(1000)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)

  # BANKSY clusters at a leiden resolution rather than a fixed k, so sweep resolutions
  # and keep those that happen to give the ground-truth number of domains.
  resolutions <- seq(0.1, 1.5, by = 0.05)
  se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = resolutions)
  cnames <- Banksy::clusterNames(se)

  cluster_counts <- sapply(cnames, function(x) {
    length(unique(colData(se)[[x]]))
  })

  valid <- cnames[cluster_counts == cluster.number]
  valid <- valid[-1]
  # Take the middle valid resolution
  if (length(valid) > 0) {
    resolution <- valid[length(valid) %/% 2 + 1]
  } else {
    resolution <- NULL
  }

  labels <- colData(se)[[resolution]]
  save(labels, file = file.path(dir.output,
                                paste0("DLPFC_BANKSY_slice_", slices[slice_index], "_labels.RData")))

  ari <- mclust::adjustedRandIndex(labels, colData(se)$clust_annotation)
  banksy_aris <- c(banksy_aris, ari)

  emb <- reducedDim(se, "PCA_M1_lam0.2")
  save(emb, file = file.path(dir.output,
                             paste0("DLPFC_BANKSY_slice_", slices[slice_index], "_embeddings.RData")))

  cat("Slice:", slices[slice_index], " resolution:", resolution, " ARI:", ari, "\n")
}

save(banksy_aris, file = file.path(dir.output, "DLPFC_BANKSY_aris.RData"))
