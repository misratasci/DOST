# BASS on the mHypothalamus (MERFISH mouse hypothalamus) dataset.
#
# Adapted from https://benchmarkst-reproducibility.readthedocs.io/en/latest/BASS_clustering.html
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
#   BASS_mHypothalamus.RData   results (per section) and aris (named list of ARIs)

library(BASS)

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"
# ---------------------------------------------------------------------------

filename <- file.path(dir.input, "MERFISH_Animal1_cnts.xlsx")
infoname <- file.path(dir.input, "MERFISH_Animal1_info.xlsx")

# Sheet IDs with domain annotations
sheets <- c('-0.04', '-0.09', '-0.14', '-0.19', '-0.24')

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

results <- list()
aris <- list()

for (sheet in sheets) {
  cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
  row.names(cnts) <- cnts[,"...1"]
  cnts <- cnts[ -c(1) ]
  cnts <- list(cnts)

  xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
  row.names(xys) <- xys[,"...1"]
  gtlabels <- xys$z
  xys <- xys[-c(1)]
  xys <- xys[-c(-2:-1)]
  xys <- list(xys)

  # C is the number of cell types, R the number of spatial domains
  C <- 20
  R <- length(unique(gtlabels))

  set.seed(1999)
  BASS <- createBASSObject(cnts, xys, C = C, R = R,
                           beta_method = "SW", init_method = "mclust",
                           nsample = 1000)

  BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
                          geneSelect = "sparkx", nSE = 3000, doPCA = TRUE,
                          scaleFeature = FALSE, nPC = 20)

  BASS <- BASS.run(BASS)
  BASS <- BASS.postprocess(BASS)

  res <- BASS@results
  results[[sheet]] <- res
  aris[[sheet]] <- mclust::adjustedRandIndex(res$z[[1]], gtlabels)
  cat("Section:", sheet, " ARI:", aris[[sheet]], "\n")
}

save(results, aris, file = file.path(dir.output, "BASS_mHypothalamus.RData"))
