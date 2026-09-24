# DOST on the mHypothalamus (MERFISH mouse hypothalamus) dataset.
#
# Two runs per section. The domain run uses R = number of annotated domains and
# lambda = 0.2, which is what the benchmark figure compares. The cell-type run
# uses R = number of annotated cell classes and lambda = 0, showing that the
# same model does cell typing when the spatial term is switched off.
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
#   mHypothalamus_DOST_sheet_<sheet>_labels.RData              labels
#   mHypothalamus_DOST_sheet_<sheet>_embeddings.RData          emb
#   mHypothalamus_DOST_aris.RData                              dost_aris, one per sheet
#   mHypothalamus_DOST_celltype_sheet_<sheet>_labels.RData     labels
#   mHypothalamus_DOST_celltype_sheet_<sheet>_embeddings.RData emb
#   mHypothalamus_DOST_celltype_aris.RData                     dost_celltype_aris

library(DOST)

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

# MERFISH settings: cells sit closer together than Visium spots, so the
# neighbourhood threshold is wider and lambda larger. The panel carries fewer
# genes than nGenes, so HVG selection keeps all of them.
neighborhood_threshold <- 3
embedding_dim <- 10
lambda <- 0.2

# Number of annotated cell classes, used by the cell-type run
n_cell_classes <- 15

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

dost_aris <- c()
dost_celltype_aris <- c()

for (sheet in sheets) {

  set.seed(1999)

  cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
  row.names(cnts) <- cnts[, "...1"]
  cnts <- cnts[-c(1)]

  xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
  row.names(xys) <- xys[, "...1"]
  gtlabels <- xys$z
  celllabels <- xys$Cell_class
  xys <- xys[-c(1)]
  xys <- xys[-c(-2:-1)]
  xys <- xys[, c(2, 1)]

  R <- length(unique(gtlabels))

  # --- spatial domains ---
  results <- DOST(cnts, xys, R,
                  neighborhood_threshold = neighborhood_threshold,
                  embedding_dim = embedding_dim,
                  lambda = lambda,
                  refinement = FALSE)

  ari <- mclust::adjustedRandIndex(results$labels, gtlabels)
  dost_aris <- c(dost_aris, ari)

  labels <- results$labels
  save(labels, file = file.path(dir.output,
                                paste0("mHypothalamus_DOST_sheet_", sheet,
                                       "_labels.RData")))
  emb <- results$Z
  save(emb, file = file.path(dir.output,
                             paste0("mHypothalamus_DOST_sheet_", sheet,
                                    "_embeddings.RData")))

  # --- cell types: as many clusters as cell classes, no spatial term ---
  results_cell <- DOST(cnts, xys, R = n_cell_classes,
                       neighborhood_threshold = neighborhood_threshold,
                       embedding_dim = embedding_dim,
                       lambda = 0,
                       refinement = FALSE)

  ari_cell <- mclust::adjustedRandIndex(results_cell$labels, celllabels)
  dost_celltype_aris <- c(dost_celltype_aris, ari_cell)

  labels <- results_cell$labels
  save(labels, file = file.path(dir.output,
                                paste0("mHypothalamus_DOST_celltype_sheet_", sheet,
                                       "_labels.RData")))
  emb <- results_cell$Z
  save(emb, file = file.path(dir.output,
                             paste0("mHypothalamus_DOST_celltype_sheet_", sheet,
                                    "_embeddings.RData")))

  cat("Section:", sheet, " domain ARI:", ari, " cell-type ARI:", ari_cell, "\n")
}

save(dost_aris, file = file.path(dir.output, "mHypothalamus_DOST_aris.RData"))
save(dost_celltype_aris,
     file = file.path(dir.output, "mHypothalamus_DOST_celltype_aris.RData"))
