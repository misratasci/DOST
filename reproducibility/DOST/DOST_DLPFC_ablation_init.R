# Ablation: what the DOST embedding is initialized with, on one DLPFC12 slice.
#
# DOST() always initializes from classical MDS of the expression distances, so
# there is no argument to vary. This script therefore runs the same pipeline
# step by step through the package internals and swaps only the starting point:
#
#   mds     classical MDS of D_expr, exactly what DOST() does
#   pca     X_ projected onto the rotation of a centered, scaled PCA of X_
#   random  i.i.d. standard normal, seeded
#
# Outputs written to dir.output:
#   DOST_DLPFC<slice_index>_init_<init>.RData   results, ari
#
# results has the same shape DOST() returns (Z, labels, losses).

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Initializations to compare
init_types <- c("random", "pca", "mds")

# Slice to run (1 to 12)
# Slice 1 and 9 were run for the ablation study
slice_index <- 1
# ---------------------------------------------------------------------------

embedding_dim <- 20

refinement <- TRUE

# Slice IDs
slices <- c(151507:151510, 151669:151676)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

cat(paste0("\nSlice ", slice_index, "\n"))

# 7 domains for slices 1-4 and 9-12, 5 for slices 5-8
if (slice_index > 4 & slice_index < 9) {
  R <- 5
} else {
  R <- 7
}

sample <- load_DLPFC_sample(slices[slice_index], dir.input)
X <- Seurat::GetAssayData(sample, layer = "counts")
coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]
N <- ncol(X)

cat("Preprocessing data...\n")
processed <- DOST:::preprocess(X, coords)
X_ <- processed$X_
D_expr <- processed$D_expr
V <- processed$V

init_random <- function(N, K, seed = 1999) {
  set.seed(seed)
  Z0 <- matrix(rnorm(N * K), nrow = N, ncol = K)
  return(Z0)
}

init_pca <- function(X, K) {
  pca_res <- prcomp(X, center = TRUE, scale. = TRUE)
  Z0 <- (X %*% pca_res$rotation)[, 1:K]
  return(Z0)
}

for (init in init_types) {

  cat("Initializing (", init, ")...\n", sep = "")
  Z_init <- switch(init,
    mds    = DOST:::mycmdscale(D_expr, embedding_dim),
    pca    = init_pca(X_, embedding_dim),
    random = init_random(N, embedding_dim),
    stop("unknown initialization: ", init)
  )

  opt <- DOST:::optimize(Z_init, D_expr, V, lambda = 0.03)
  Z <- opt$Z
  losses <- opt$losses

  max_mclust_c <- DOST:::cluster_embedding(Z, R, modelName = "EEE")
  if (is.null(max_mclust_c)) {
    cat("Mclust failed!\n")
    next
  }

  if (refinement) {
    cat("Refining labels...\n")
    labels <- DOST:::refine_labels(coords, max_mclust_c$classification)
  } else {
    labels <- max_mclust_c$classification
  }

  results <- list(Z = Z, labels = labels, losses = losses)
  ari <- mclust::adjustedRandIndex(labels, sample@meta.data$layers)

  cat("slice:", slices[slice_index], " init:", init, " ARI:", ari, "\n")

  save(results, ari,
       file = file.path(dir.output,
                        paste0("DOST_DLPFC", slice_index, "_init_", init, ".RData")))
}
