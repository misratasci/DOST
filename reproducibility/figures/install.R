
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c("ggplot2", "patchwork", "scales", "reshape2", "ggbeeswarm",
                   "dplyr", "tidyr", "knitr", "umap", "mclust", "readxl",
                   "Seurat", "cowplot", "kableExtra"))

BiocManager::install("hdf5r")
