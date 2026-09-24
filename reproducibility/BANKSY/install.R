
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# CRAN
install.packages(c("Seurat", "SeuratObject", "mclust", "hdf5r", "readxl"))

# Bioconductor
BiocManager::install(c("Banksy", "SpatialExperiment", "SummarizedExperiment",
                       "scuttle", "scater"))
