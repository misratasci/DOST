
cluster_embedding <- function(Z, R, modelName = "EEE") {
  cat("\nClustering...\n")
  mclust_c <- mclust::Mclust(data = Z, G = R, modelNames = modelName, verbose = FALSE, subset = Inf)
  if (is.null(mclust_c)) {
    cat("Mclust failed, retrying with different model families...\n")
    mclust_c <- mclust::Mclust(data = Z, G = R, verbose = FALSE, subset = Inf)
  }
  if (is.null(mclust_c)) {
    cat("Mclust failed\n")
  }
  return (mclust_c)
}

refine_labels <- function(coords, labels, k = 30) {
  n <- nrow(coords)
  if (k >= n) stop("k must be < number of rows in Z")

  D <- compute_D(coords)

  new_labels <- character(n)
  for (i in 1:n) {
    ord <- order(D[i, ], decreasing = FALSE)
    neigh <- ord[ord != i][1:k]
    tab <- table(labels[neigh])
    new_labels[i] <- names(which.max(tab))
  }
  return(new_labels)
}
