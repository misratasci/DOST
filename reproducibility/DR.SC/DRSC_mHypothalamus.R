# DR.SC on the mHypothalamus (MERFISH mouse hypothalamus) dataset.
#
# Method: https://feiyoung.github.io/DR.SC/
# Data loading from BenchmarkST
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
#   mHypothalamus_DRSC_sheet_<sheet>_labels.RData       labels, the chosen clustering
#   mHypothalamus_DRSC_sheet_<sheet>_embeddings.RData   emb, the DR.SC embedding
#   mHypothalamus_DRSC_aris.RData                drsc_aris, one ARI per sheet in the loop

library(Seurat)
library(DR.SC)

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir.output <- "path/to/output"
# ---------------------------------------------------------------------------

filename <- file.path(dir.input, "MERFISH_Animal1_cnts.xlsx")
infoname <- file.path(dir.input, "MERFISH_Animal1_info.xlsx")

# Sheet IDs with domain annotations
sheets <- c('-0.04', '-0.09', '-0.14', '-0.19', '-0.24')

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

drsc_aris <- c()

for (sheet in sheets) {
  cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
  row.names(cnts) <- cnts[,"...1"]
  cnts <- cnts[ -c(1) ]

  xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
  row.names(xys) <- xys[,"...1"]
  gtlabels <- xys$z
  xys <- xys[-c(1)]

  sample <- CreateSeuratObject(counts = cnts, project = "43F", min.cells = 3, names.delim = "-", names.field = 2)

  sample <- AddMetaData(sample,
                         metadata = xys$x,
                         col.name = 'row')
  sample <- AddMetaData(sample,
                         metadata = xys$y,
                         col.name = 'col')
  sample <- AddMetaData(sample,
                         metadata = xys$z,
                         col.name = 'layer_guess_reordered')

  sample$orig.ident <- 1
  Idents(sample) <- row.names(sample@meta.data)

  cluster.number <- length(unique(gtlabels))

  # No gene selection for MERFISH
  sample <- NormalizeData(sample, verbose = FALSE)
  seu <- FindVariableFeatures(sample, nfeatures = 500, verbose = F)

  seu <- DR.SC(seu, K = as.numeric(cluster.number), platform = 'Other_SRT', verbose = FALSE)

  ari_drsc <- mclust::adjustedRandIndex(seu$spatial.drsc.cluster, gtlabels)
  print(paste0("ARI: ", ari_drsc))
  drsc_aris <- c(drsc_aris, ari_drsc)

  emb <- Embeddings(seu, reduction = "dr-sc")
  save(emb, file = file.path(dir.output,
                             paste0("mHypothalamus_DRSC_sheet_", sheet, "_embeddings.RData")))
  labels <- seu$spatial.drsc.cluster
  save(labels, file = file.path(dir.output,
                                paste0("mHypothalamus_DRSC_sheet_", sheet, "_labels.RData")))

}

save(drsc_aris, file = file.path(dir.output, "mHypothalamus_DRSC_aris.RData"))
