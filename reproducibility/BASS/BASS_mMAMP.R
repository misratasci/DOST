# BASS on the mMAMP (mouse brain anterior, "MA") dataset.
#
# Adapted from https://benchmarkst-reproducibility.readthedocs.io/en/latest/BASS_clustering.html
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
#   mMAMP_BASS_labels.RData   zlabels (domains), clabels (cell types) and the ARI

library(Seurat)
library(BASS)
source("load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mMAMP"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"
# ---------------------------------------------------------------------------

section <- "MA"

set.seed(1999)
dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

sample <- load_mMAMP_sample(dir.input, section = section)
cnts <- list(sample@assays$Spatial$counts)

df_meta <- read.table(file.path(dir.input, section, "gt", "tissue_positions_list_GTs.txt"),
                      sep = "\t", header = TRUE)
xym <- data.frame(row = df_meta$array_row, col = df_meta$array_col)
row.names(xym) <- row.names(sample@meta.data)
xym <- list(xym)

C <- length(unique(sample@meta.data$ground_truth))
R <- length(unique(sample@meta.data$ground_truth))

BASS <- createBASSObject(cnts, xym, C = C, R = R)
listAllHyper(BASS)
BASS <- BASS.preprocess(BASS,
                        geneSelect = "sparkx" # or "hvgs"
)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS)

zlabels <- BASS@results$z
clabels <- BASS@results$c

gt <- sample@meta.data$ground_truth
ari <- mclust::adjustedRandIndex(zlabels[[1]], gt)
cat("Dataset: mMAMP", section, " ARI:", ari, "\n")

save(zlabels, clabels, ari, file = file.path(dir.output, "mMAMP_BASS_labels.RData"))
