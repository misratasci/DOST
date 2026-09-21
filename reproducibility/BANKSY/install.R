# Package installation for the BANKSY reproducibility scripts.
#
#   Rscript install.R
#
# Banksy is on Bioconductor (as "Banksy"). The scripts also use the SpatialExperiment
# / SummarizedExperiment stack and Seurat for normalization and feature selection.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# CRAN
install.packages(c("Seurat", "SeuratObject", "mclust", "hdf5r"))

# Bioconductor
BiocManager::install(c("Banksy", "SpatialExperiment", "SummarizedExperiment",
                       "scuttle", "scater"))
