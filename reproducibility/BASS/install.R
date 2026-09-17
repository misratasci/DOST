if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# CRAN
install.packages(c("Seurat", "SeuratObject", "mclust", "readxl", "hdf5r"))

# GitHub
remotes::install_github("zhengli09/BASS")
remotes::install_github("xzhoulab/SPARK")
