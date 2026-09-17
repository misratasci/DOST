# BASS on the BC (human breast cancer, Visium) dataset.
#
# Adapted from https://benchmarkst-reproducibility.readthedocs.io/en/latest/BASS_clustering.html
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download BC data from https://zenodo.org/records/10698903
#
# Expected layout:
#   BC/section1/{spatial/, section1_filtered_feature_bc_matrix.h5,
#                gt/gold_metadata.tsv, gt/tissue_positions_list_GTs.txt}
#
# Outputs written to dir.output:
#   BC_BASS_labels.RData   zlabels (domains), clabels (cell types) and the ARI

library(Seurat)
library(BASS)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/BC"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"
# ---------------------------------------------------------------------------

set.seed(1999)
dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_BC_sample(dir.input)
cnts <- list(sample@assays$Spatial$counts)

df_meta <- read.table(file.path(dir.input, "section1", "gt", "tissue_positions_list_GTs.txt"),
                      sep = ",", row.names = 1)
xym <- data.frame(row = df_meta['V3'], col = df_meta['V4'])
row.names(xym) <- row.names(sample@meta.data)
xym <- list(xym)

# C is the number of cell types, R the number of spatial domains
C <- 50
R <- 20

BASS <- createBASSObject(cnts, xym, C = C, R = R)
listAllHyper(BASS)
BASS <- BASS.preprocess(BASS,
                        geneSelect = "sparkx" # or "hvgs"
)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS)

zlabels <- BASS@results$z
clabels <- BASS@results$c

gt <- sample@meta.data$fine_annot_type
ari <- mclust::adjustedRandIndex(zlabels[[1]], gt)
cat("Dataset: BC section1  ARI:", ari, "\n")

save(zlabels, clabels, ari, file = file.path(dir.output, "BC_BASS_labels.RData"))
