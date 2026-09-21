# BANKSY on the mHypothalamus (MERFISH mouse hypothalamus) dataset.
#
# Method: https://prabhakarlab.github.io/Banksy/
# Preprocessing follows https://prabhakarlab.github.io/Banksy/articles/multi-sample.html
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download mHypothalamus data from https://zenodo.org/records/10698909
#
# Expected layout:
#   mHypothalamus/MERFISH_Animal1_cnts.xlsx   (one sheet per section)
#   mHypothalamus/MERFISH_Animal1_info.xlsx   (same sheet names)
#
# Outputs written to dir.output:
#   mHypothalamus_BANKSY_sheet_<sheet>_labels.RData       labels, the chosen clustering
#   mHypothalamus_BANKSY_sheet_<sheet>_embeddings.RData   emb, the BANKSY PCA embedding
#   mHypothalamus_BANKSY_aris.RData                banksy_aris, one ARI per sheet in the loop

library(Banksy)

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "/Users/misratasci/Desktop/ST-research/data/mHypothalamus"  #"path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir.output <- "output" #"path/to/output/"
# ---------------------------------------------------------------------------

filename <- file.path(dir.input, "MERFISH_Animal1_cnts.xlsx")
infoname <- file.path(dir.input, "MERFISH_Animal1_info.xlsx")

# Sheet IDs with domain annotations
sheets <- c('-0.04', '-0.09', '-0.14', '-0.19', '-0.24')

# BANKSY default hyperparameters
lambda <- c(0.8)
k_geom <- c(15, 30)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

results <- list()
banksy_aris <- c()

for (sheet in sheets) {
  gcm <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
  row.names(gcm) <- gcm[,"...1"]
  gcm <- gcm[ -c(1) ]
  gcm <- as.matrix(gcm)

  locs <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
  row.names(locs) <- locs[,"...1"]
  gtlabels <- locs$z
  locs <- locs[-c(1)]
  locs <- locs[-c(-2:-1)]
  locs <- as.matrix(locs)

  cluster.number <- length(unique(gtlabels))

  se <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)

  colData(se) <- DataFrame(
    sample_id = sheet,
    clust_annotation = factor(
      addNA(gtlabels),
      exclude = NULL, labels = seq(cluster.number)
    ),
    row.names = row.names(locs)
  )

  # No HVG selection for MERFISH, just normalize and run BANKSY
  seu <- as.Seurat(se, data = NULL)
  seu <- NormalizeData(seu, scale.factor = 5000, normalization.method = "RC")
  aname <- "normcounts"
  assay(se, aname) <- GetAssayData(seu)

  rm(seu)

  se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

  set.seed(1000)
  se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
  se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)

  # BANKSY clusters at a leiden resolution rather than a fixed k, so sweep resolutions
  # and keep those that happen to give the ground-truth number of domains.
  resolutions <- seq(0.1, 2.5, by = 0.01)
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
                                paste0("mHypothalamus_BANKSY_sheet_", sheet, "_labels.RData")))

  ari <- mclust::adjustedRandIndex(labels, colData(se)$clust_annotation)
  banksy_aris <- c(banksy_aris, ari)

  emb <- reducedDim(se, "PCA_M1_lam0.8")
  save(emb, file = file.path(dir.output,
                             paste0("mHypothalamus_BANKSY_sheet_", sheet, "_embeddings.RData")))

  cat("Slice:", sheet, " resolution:", resolution, " ARI:", ari, "\n")

}

save(results, banksy_aris, file = file.path(dir.output, "mHypothalamus_BANKSY_aris.RData"))
