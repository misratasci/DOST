
get_HVG <- function(X, nHVG = 3000) {
  dec <- scran::modelGeneVar(X)
  hvg <- scran::getTopHVGs(dec, n = nHVG)
  gene_idx <- match(hvg, rownames(X))
  gene_idx <- gene_idx[!is.na(gene_idx)]
  X_hvg <- X[gene_idx, ]
  #sparse matrix (genes x spots) to dense (spots x genes)
  X_hvg <- t(as.matrix(X_hvg))
  return (X_hvg)
}

get_SVG <- function(X, coords, nSVG = 3000) {
  if (!requireNamespace("SPARK", quietly = TRUE)) {
    stop(
      "The 'SPARK' package is required for SVG selection. ",
      "Please install it using: devtools::install_github('xzhoulab/SPARK')",
      call. = FALSE
    )
  }
  sparkx <- SPARK::sparkx(X, coords, verbose = FALSE)
  sparkx <- sparkx$res_mtest[order(sparkx$res_mtest$adjustedPval), ]
  svg <- head(rownames(sparkx), n = nSVG)
  gene_idx <- match(svg, rownames(X))
  gene_idx <- gene_idx[!is.na(gene_idx)]
  X_svg <- X[gene_idx, ]
  #sparse matrix (genes x spots) to dense (spots x genes)
  X_svg <- t(as.matrix(X_svg))
  return (X_svg)
}

build_adj_mat <- function(coords, threshold_level = 1) {
  N <- nrow(coords)
  V <- matrix(0, nrow = N, ncol = N)
  dist <- as.matrix(dist(coords, method = "euclidean"))
  dist_no_self <- dist
  diag(dist_no_self) <- Inf
  r1 <- apply(dist_no_self, 1, min)
  threshold_level <- 1.5 * median(r1) * threshold_level
  V[dist <= threshold_level] <- 1
  return (V)
}

preprocess <- function(X, coords, selected_genes = "HVG", nGenes = 3000, neighborhood_threshold = 1) {
  #log normalize
  X <- scuttle::normalizeCounts(X, log = TRUE)
  #each cell (spot) is divided by the total read count of that cell, then log2 transform

  #get highly variable and spatially variable genes
  if (selected_genes == "HVG") {
    cat("Selecting HVGs\n")
    X_ <- get_HVG(X, nGenes)
  } else if (selected_genes == "SVG") {
    cat("Selecting SVGs\n")
    X_ <- get_SVG(X, coords, nGenes)
  } else {
    cat("Using all genes\n")
    X_ <- t(as.matrix(X))
  }

  D_expr <- compute_D(X_)
  K_expr <- exp(- D_expr^2 / (0.5 * mean(D_expr)^2))
  D_expr <- sqrt(2 - 2 * K_expr)

  #build neighborhood graph V
  V <- build_adj_mat(coords, threshold_level = neighborhood_threshold)

  return (list(X_ = X_, D_expr = D_expr, V = V))
}

compute_D <- function(X) {
  D <- Rfast::Dist(X, method = "euclidean", result = "matrix")
  return (D)
}

