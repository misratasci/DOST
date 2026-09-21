# BANKSY on the mMAMP (mouse brain anterior, "MA") dataset.
#
# Method: https://prabhakarlab.github.io/Banksy/
# Preprocessing follows https://prabhakarlab.github.io/Banksy/articles/multi-sample.html
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download mMAMP data from https://zenodo.org/records/10698931
#
# Expected layout:
#   mMAMP/MA/{spatial/, MA_filtered_feature_bc_matrix.h5, metadata.tsv,
#             gt/tissue_positions_list_GTs.txt}
#
# Outputs written to dir.output:
#   mMAMP_BANKSY_labels.RData              labels, the chosen clustering
#   mMAMP_BANKSY_embeddings.RData          emb, the BANKSY embedding

library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(Seurat)
library(hdf5r)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mMAMP"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"
# ---------------------------------------------------------------------------

# BANKSY recommends the following hyperparameters for 10x Visium v1v2
# https://prabhakarlab.github.io/Banksy/articles/parameter-selection.html
lambda <- c(0.2)
k_geom <- c(18, 18)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_mMAMP_sample(dir.input)
cluster.number <- length(unique(sample@meta.data$ground_truth))
gcm <- LayerData(sample, assay = "Spatial", layer = "counts")

locs <- data.frame(sdimx = Seurat::GetTissueCoordinates(sample)$x, sdimy = Seurat::GetTissueCoordinates(sample)$y)
row.names(locs) <- row.names(sample@meta.data)
spatial_coor <- as.matrix(locs)
rownames(spatial_coor) <- row.names(locs)

se <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = spatial_coor)

colData(se) <- DataFrame(
  sample_id = "mMAMP",
  clust_annotation = factor(
    addNA(sample@meta.data$ground_truth),
    exclude = NULL, labels = seq(length(unique(sample@meta.data$ground_truth)))
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
resolutions <- seq(4.5, 10.5, by = 0.05)
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
                              paste0("mMAMP_BANKSY_labels.RData")))

ari <- mclust::adjustedRandIndex(labels, colData(se)$clust_annotation)
print(paste0("ARI: ", ari))
emb <- reducedDim(se, "PCA_M1_lam0.2")
save(emb, file = file.path(dir.output,
                           paste0("mMAMP_BANKSY_embeddings.RData")))
