# Installs what the DOST reproducibility scripts need.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# DOST's own dependencies, plus what the scripts use to load the data
BiocManager::install(c("scran", "scuttle"))
install.packages(c("Rcpp", "RSpectra", "Rfast", "mclust", "Seurat", "readxl"))

# hdf5r is what Seurat::Load10X_Spatial needs to read the Visium .h5 matrices
BiocManager::install("hdf5r")

# SPARK is only needed for the SVG arm of DOST_DLPFC_ablation_nGenes.R
devtools::install_github("xzhoulab/SPARK")

# DOST itself, from the root of this repository
devtools::install_local("../..")
