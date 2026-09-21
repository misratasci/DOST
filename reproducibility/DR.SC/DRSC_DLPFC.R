# DR.SC on the DLPFC12 dataset.
#
# Method: https://feiyoung.github.io/DR.SC/
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download DLPFC12 data from https://zenodo.org/records/10698880
#
# Outputs written to dir.output:
#   DLPFC_DRSC_labels_<slice>.RData       labels, the spatial.drsc.cluster assignment
#   DLPFC_DRSC_embeddings_<slice>.RData   emb, the "dr-sc" reduction
#   DLPFC_DRSC_aris.RData                 drsc_aris, one ARI per slice in the loop

library(Seurat)
library(DR.SC)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"

# Slice indices to run (1 to 12).
slice_indices <- 1:12
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

# Number of highly variable genes DR.SC is given
n_features <- 500

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

drsc_aris <- c()

for (slice_index in slice_indices) {

  if (slice_index > 4 & slice_index < 9) {
    cluster.number <- 5
  } else {
    cluster.number <- 7
  }

  sample <- load_DLPFC_sample(slices[slice_index], dir.input)
  sample <- NormalizeData(sample, verbose = FALSE)
  seu <- FindVariableFeatures(sample, nfeatures = n_features, verbose = FALSE)

  seu <- DR.SC(seu, K = as.numeric(cluster.number), platform = 'Visium', verbose = FALSE)

  ari_drsc <- mclust::adjustedRandIndex(seu$spatial.drsc.cluster, seu@meta.data$layers)
  drsc_aris <- c(drsc_aris, ari_drsc)

  emb <- Embeddings(seu, reduction = "dr-sc")
  save(emb, file = file.path(dir.output,
                             paste0("DLPFC_DRSC_embeddings_", slices[slice_index], ".RData")))
  labels <- seu$spatial.drsc.cluster
  save(labels, file = file.path(dir.output,
                                paste0("DLPFC_DRSC_labels_", slices[slice_index], ".RData")))

  cat("Slice:", slices[slice_index], " ARI:", ari_drsc, "\n")
}

save(drsc_aris, file = file.path(dir.output, "DLPFC_DRSC_aris.RData"))
