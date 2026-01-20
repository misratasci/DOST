#' @importFrom Rcpp evalCpp
#' @import mclust
#' @useDynLib DOST, .registration = TRUE
NULL

mycmdscale <- function (d, k) {
  if (is.null(n <- attr(d, "Size"))) {
    x <- as.matrix(d^2)
    storage.mode(x) <- "double"
    if ((n <- nrow(x)) != ncol(x))
      stop("distances must be result of 'dist' or a square matrix")
    rn <- rownames(x)
  } else {
    rn <- attr(d, "Labels")
    x <- matrix(0, n, n)
    x[row(x) > col(x)] <- d^2
    x <- x + t(x)
  }
  n <- as.integer(n)
  if (is.na(n) || n > 46340)
    stop(gettextf("invalid value of %s", "'n'"), domain = NA)
  if ((k <- as.integer(k)) > n - 1 || k < 1)
    stop("'k' must be in {1, 2, ..  n - 1}")
  R = x*0 + rowMeans(x)
  C = t(x*0 + colMeans(x))
  x = x - R - C + mean(x[])
  e <- RSpectra::eigs_sym(-x/2, k = k, which = "LM")
  ev <- e$values
  evec <- e$vectors[, seq_len(k), drop = FALSE]
  k1 <- sum(ev > 0)
  if (k1 < k) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0",
                     k1, k), domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
  }
  points <- evec * rep(sqrt(ev), each = n)
  dimnames(points) <- list(rn, NULL)
  return (points)
}

compute_loss <- function(D_lab, D_expr, V, lambda, eps = 1e-20) {
  N <- nrow(D_lab)
  L1 <- sum( (D_lab - D_expr)^2 / (D_expr + diag(eps, N, N)) ) / nrow(D_lab)^2
  L2 <- sum(V * D_lab) / sum(V)
  return (L1 + lambda * L2)
}

optimize <- function(Z, D_expr, V, lambda, lr = 16, eps = 1e-20, max_iterations = 20, loss_tol = 1e-5) {
  D_lab <- compute_D(Z)
  current_loss <- compute_loss(D_lab, D_expr, V, lambda, eps)
  losses <- c(current_loss)
  i <- 0
  N <- nrow(D_lab)
  grad_a <- 2 / (N^2 * D_expr)
  grad_b <- ((lambda * V) / sum(V)) - (2 / N^2)
  Zs <- list()
  while (i < max_iterations) {
    cat("\rOptimizing - Iteration:", i + 1, "of" , max_iterations, " Loss:", format(current_loss, digits = 6))
    grad <- compute_grad(Z, D_lab, grad_a, grad_b, eps)
    coef <- 1
    best_loss <- current_loss
    best_Z <- Z
    while (TRUE) {
      new_Z <- Z - coef * lr * grad
      new_D_lab <- compute_D(new_Z)
      modified_loss <- compute_loss(new_D_lab, D_expr, V, lambda, eps)
      if (modified_loss <= best_loss) {
        if (log2(coef) < 0) {
          best_loss <- modified_loss
          Z <- best_Z
          break
        }
        else {
          coef <- coef * 2
          best_loss <- modified_loss
          best_Z <- new_Z
        }
      }
      else {
        if (log2(coef) <= 0 && coef > 1e-3) {
          coef <- coef / 2
        }
        else {
          coef <- coef / 2
          Z <- best_Z
          break
        }
      }
    }
    D_lab <- compute_D(Z)
    new_loss <- compute_loss(D_lab, D_expr, V, lambda, eps)
    losses <- c(losses, new_loss)
    if ((current_loss - new_loss)/current_loss < loss_tol) {
      print("Loss improvement is too small, stopping optimization.")
      break
    }
    current_loss <- new_loss
    i <- i + 1
    Zs[[i]] <- Z
  }
  return(list(Zs = Zs, Z = Z, losses = losses))
}

#' DOST: Distance-preserving Optimization for Spatial Transcriptomics
#'
#' Runs DOST to identify spatial domains in spatial transcriptomics data.
#'
#' @param X A sparse matrix or matrix of gene expression counts (Genes x Spots).
#' @param coords A data frame or matrix of spatial coordinates (Spots x 2).
#'   Rows of \code{coords} must match the columns of \code{X}.
#' @param R Integer. The number of clusters (domains) to identify.
#' @param selected_genes Character. Gene selection method: \code{"HVG"} (default), \code{"SVG"}, or \code{"all"}.
#' @param nGenes Integer. Number of genes to use if selection is enabled. Default is 3000.
#' @param refinement Logical. Whether to smooth the final labels using spatial neighbors.
#'   Recommended for tissues with laminar structures (e.g., DLPFC, mPFC). Default is \code{TRUE}.
#'
#' @param neighborhood_threshold Spatial adjacency threshold multiplier. Default is 1.
#' @param embedding_dim Dimension of the initial embedding. Default is 20.
#' @param lambda Regularization weight for spatial coherence. Default is 0.03.
#' @param lr Optimization learning rate. Default is 16.
#' @param max_iterations Max optimization steps. Default is 20.
#' @param refine_k Number of neighbors for label refinement. Default is 30.
#'
#' @details
#' \strong{Recommended Settings by Technology:}
#'
#' The default parameters are optimized for **10x Visium** data. For other technologies,
#' we recommend the following adjustments based on our benchmarks:
#'
#' \itemize{
#'   \item \strong{10x Visium:}
#'     \itemize{
#'       \item \code{embedding_dim = 20}
#'       \item \code{neighborhood_threshold = 1}
#'       \item \code{lambda = 0.03}
#'       \item \code{nGenes = 3000} (HVG or SVG)
#'     }
#'   \item \strong{MERFISH:}
#'     \itemize{
#'       \item \code{embedding_dim = 10}
#'       \item \code{neighborhood_threshold = 3} (Due to higher spatial resolution)
#'       \item \code{lambda = 0.2}
#'       \item \code{selected_genes = "all"} (If gene count is low)
#'     }
#'   \item \strong{STARmap:}
#'     \itemize{
#'       \item \code{embedding_dim = 10}
#'       \item \code{neighborhood_threshold = 1}
#'       \item \code{lambda = 0.15}
#'       \item \code{selected_genes = "all"}
#'     }
#' }
#'
#' \strong{Refinement Note:}
#' The optional refinement step (\code{refinement = TRUE}) effectively smooths boundaries
#' and is highly recommended for tissues with known laminar or layered structures
#' (e.g., Cortex). For tissues with discrete, punctate domains, you may wish to
#' set this to \code{FALSE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Z}: The low-dimensional embedding matrix.
#'   \item \code{labels}: The final cluster labels.
#'   \item \code{losses}: The loss function history.
#' }
#' @export
DOST <- function(X, coords, R, selected_genes = "HVG",
                 nGenes = 3000, neighborhood_threshold = 1, embedding_dim = 20,
                 lambda = 0.03, lr = 16, max_iterations = 20,
                 refinement = TRUE, refine_k = 30) {
  cat("Preprocessing data...\n")
  processed <- preprocess(X, coords, selected_genes, nGenes, neighborhood_threshold)
  D_expr <- processed$D_expr
  V <- processed$V
  cat("Initializing...\n")
  Z_init <- mycmdscale(D_expr, embedding_dim)
  results <- optimize(Z_init, D_expr, V, lambda, lr = lr, eps = 1e-20,
                      max_iterations = max_iterations, loss_tol = 1e-5)
  Z <- results$Z
  losses <- results$losses
  max_mclust_c <- cluster_embedding(Z, R, modelName = "EEE")
  if (is.null(max_mclust_c)) {
    cat("Mclust failed!\n")
    return (NULL)
  }
  if (refinement) {
    cat("Refining labels...\n")
    labels <- refine_labels(coords, max_mclust_c$classification, refine_k)
  } else {
    labels <- max_mclust_c$classification
  }
  return (list(Z = Z, losses = losses, labels = labels))
}
