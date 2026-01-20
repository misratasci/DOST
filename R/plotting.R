
utils::globalVariables(c("x", "y", "cluster", "UMAP1", "UMAP2"))

#' Plot DOST Spatial Clusters
#'
#' Visualizes the clustering results in spatial coordinates.
#'
#' @param coords A data frame or matrix of spatial coordinates (Spots x 2).
#' @param labels A vector of cluster labels.
#' @param point_size Numeric, size of the points (default 1).
#' @return A ggplot object.
#' @export
plot_DOST_spatial <- function(coords, labels, point_size = 1) {
  df <- as.data.frame(coords)
  colnames(df)[1:2] <- c("x", "y")
  df$cluster <- as.factor(labels)

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::theme_void() +
    ggplot2::labs(color = "Cluster") +
    ggplot2::scale_color_discrete()
}

#' Plot DOST Embedding (UMAP)
#'
#' Projects the DOST embedding into 2D using UMAP.
#'
#' @param Z The embedding matrix from DOST.
#' @param labels A vector of cluster labels.
#' @export
plot_DOST_umap <- function(Z, labels) {
  umap_res <- uwot::umap(Z, n_neighbors = 15, min_dist = 0.1)
  df <- data.frame(UMAP1 = umap_res[,1], UMAP2 = umap_res[,2], cluster = as.factor(labels))

  ggplot2::ggplot(df, ggplot2::aes(x = UMAP1, y = UMAP2, color = cluster)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::theme_void()
}
