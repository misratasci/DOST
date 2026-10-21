if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# CRAN
install.packages(c("Seurat", "SeuratObject", "mclust", "readxl", "hdf5r"))

# GitHub
# BASS requires an older version of RcppArmadillo.
# To avoid issues, we install this specific version.
remotes::install_version(
  "RcppArmadillo",
  version = "14.6.3-1",
  repos = "https://cloud.r-project.org",
  upgrade = "never"
)

remotes::install_github("zhengli09/BASS")
remotes::install_github("xzhoulab/SPARK")
