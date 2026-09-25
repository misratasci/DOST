
# Helpers shared by the figure scripts.
#
# Two groups:
#   * readers   turn the files the method scripts write into label vectors,
#               embeddings and ARI vectors
#   * plotters  the panel styles the article uses

library(ggplot2)
library(patchwork)
library(scales)
library(Seurat)
library(umap)

# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

align_to_ids <- function(values, ids, target_ids, what = "labels") {
  if (!is.null(ids) && !anyDuplicated(ids) && all(target_ids %in% ids)) {
    return(values[match(target_ids, ids)])
  }
  if (length(values) != length(target_ids)) {
    stop(sprintf("%s has %d entries but %d were expected, and the ids do not match",
                 what, length(values), length(target_ids)))
  }
  warning(sprintf("could not match %s on ids, assuming file order matches the object",
                  what), call. = FALSE)
  values
}

read_labels_csv <- function(path, column, target_ids = NULL) {
  df <- utils::read.csv(path, row.names = 1)
  values <- df[[column]]
  if (is.null(target_ids)) return(values)
  align_to_ids(values, rownames(df), target_ids, basename(path))
}

read_labels_txt <- function(path) {
  utils::read.csv(path, header = FALSE)[[1]]
}

read_embedding_csv <- function(path, target_ids = NULL, ids = NULL) {
  emb <- as.matrix(utils::read.csv(path, header = FALSE))
  if (is.null(target_ids) || is.null(ids)) return(emb)
  if (nrow(emb) != length(ids)) {
    stop(sprintf("%s has %d rows but %d ids were given",
                 basename(path), nrow(emb), length(ids)))
  }
  idx <- match(target_ids, ids)
  if (anyNA(idx)) {
    stop(sprintf("%s is missing %d of the spots in the object",
                 basename(path), sum(is.na(idx))))
  }
  emb[idx, , drop = FALSE]
}

load_object <- function(path, name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!name %in% ls(env)) {
    stop(sprintf("'%s' is not in %s (it holds: %s)",
                 name, basename(path), paste(ls(env), collapse = ", ")))
  }
  get(name, envir = env)
}

# ---------------------------------------------------------------------------
# Plotters
# ---------------------------------------------------------------------------

make_spatial_plot <- function(sample, title, column, gt_column,
                              show_ari = TRUE, cols = NULL, pt.size.factor = 2.5) {
  p <- Seurat::SpatialDimPlot(sample, group.by = column,
                              pt.size.factor = pt.size.factor, alpha = 1,
                              cols = cols)
  if (show_ari) {
    ari <- mclust::adjustedRandIndex(sample@meta.data[[column]],
                                     sample@meta.data[[gt_column]])
    p <- p + labs(caption = paste0("ARI=", formatC(ari, digits = 2, format = "f")))
  }
  p +
    Seurat::NoLegend() +
    theme(plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0.5),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(title)
}

make_spatial_column <- function(samples, title, column, gt_column,
                                show_ari = TRUE, cols = NULL,
                                captions = NULL, pt.size.factor = 2.5) {
  plots <- list()
  for (i in seq_along(samples)) {
    p <- Seurat::SpatialDimPlot(samples[[i]], group.by = column,
                                pt.size.factor = pt.size.factor, alpha = 1,
                                cols = cols)
    if (!is.null(captions)) {
      p <- p + labs(caption = captions[i])
    } else if (show_ari) {
      ari <- mclust::adjustedRandIndex(samples[[i]]@meta.data[[gt_column]],
                                       samples[[i]]@meta.data[[column]])
      p <- p + labs(caption = paste0("ARI=", formatC(ari, digits = 2, format = "f")))
    }
    p <- p +
      Seurat::NoLegend() +
      theme(plot.caption.position = "plot",
            plot.caption = element_text(hjust = 0.5),
            plot.title.position = "plot",
            plot.title = element_text(hjust = 0.5))
    if (i == 1) p <- p + ggtitle(title)
    plots[[i]] <- p
  }
  wrap_plots(plots, ncol = 1)
}

make_umap_plot <- function(emb, gtlabels, title, legend = FALSE,
                           trim_quant = 0, colors = c(0, 360)) {
  sample.umap <- umap::umap(emb)
  layout <- as.data.frame(sample.umap$layout)
  colnames(layout) <- c("x", "y")

  layer_levels <- sort(unique(gtlabels))
  df <- transform(layout, layer = factor(gtlabels, levels = layer_levels))

  if (trim_quant > 0) {
    x_q <- quantile(df$x, probs = c(trim_quant, 1 - trim_quant))
    y_q <- quantile(df$y, probs = c(trim_quant, 1 - trim_quant))
    df <- df[df$x >= x_q[1] & df$x <= x_q[2] &
               df$y >= y_q[1] & df$y <= y_q[2], , drop = FALSE]
  }

  num_map <- stats::setNames(as.character(seq_along(layer_levels)), layer_levels)
  df$layer_num <- factor(num_map[as.character(df$layer)],
                         levels = as.character(seq_along(layer_levels)))
  cols <- hue_pal(h = colors)(length(layer_levels))

  p <- ggplot(df, aes(x = x, y = y, color = layer_num)) +
    geom_point(size = 0.6) +
    scale_color_manual(values = cols, labels = layer_levels) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    theme_void(base_size = 10) +
    theme(plot.margin = unit(c(0.2, 0.7, 1.2, 0.7), "lines"),
          legend.position = "right",
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(title)
  if (!legend) p <- p + Seurat::NoLegend()
  p
}

make_points_plot <- function(coords, labels, title = NULL, caption = NULL,
                             size = 0.1, legend = FALSE, legend_title = NULL,
                             colors = c(0, 360), legend_ncol = 2) {
  df <- data.frame(x = coords[, 1], y = coords[, 2], Domain = labels)
  cols <- hue_pal(h = colors)(length(unique(labels)))
  p <- ggplot(df, aes(x, y, color = Domain)) +
    geom_point(size = size) +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0.5),
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  if (!legend) p <- p + theme(legend.position = "none")
  if (legend) p <- p +
    guides(color = guide_legend(title = legend_title, ncol = legend_ncol)) +
    theme(legend.position = "bottom")
  if (!is.null(caption)) p <- p + labs(caption = caption)
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

make_points_column <- function(coords_list, labels_list, gtlabels_list, sheets,
                               title, size = 0.1, show_ari = TRUE,
                               captions = NULL) {
  plots <- list()
  for (i in seq_along(sheets)) {
    caption <- NULL
    if (!is.null(captions)) {
      caption <- captions[i]
    } else if (show_ari) {
      ari <- mclust::adjustedRandIndex(labels_list[[i]], gtlabels_list[[i]])
      caption <- paste0("ARI=", formatC(ari, digits = 2, format = "f"))
    }
    plots[[i]] <- make_points_plot(coords_list[[i]], labels_list[[i]],
                                   title = if (i == 1) title else NULL,
                                   caption = caption,
                                   size = size)
  }
  plots
}

make_ari_violin <- function(aris, show_mean = TRUE) {
  long <- reshape2::melt(aris, variable.name = "method", value.name = "ARI")
  p <- ggplot(long, aes(x = method, y = ARI, fill = method)) +
    geom_violin(trim = FALSE, color = "black", width = 0.8, alpha = 1, adjust = 1, scale = "width") +
    ggbeeswarm::geom_beeswarm(cex = 1.0, size = 1.5, priority = "density",
                              alpha = 1, color = "black")
  if (show_mean) {
    p <- p + stat_summary(fun = mean, geom = "crossbar",
                          aes(ymin = after_stat(y), ymax = after_stat(y)),
                          width = 0.6, color = "darkred", linewidth = 0.4)
  }
  p +
    guides(fill = "none", color = "none") +
    labs(y = "ARI", x = NULL) +
    theme_classic()
}

save_cropped <- function(path, plot, width, height) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, plot, width = width, height = height)
  if (nzchar(Sys.which("pdfcrop"))) {
    system2("pdfcrop", args = c(path, path))
  } else {
    message("pdfcrop not found, leaving ", basename(path), " uncropped")
  }
  invisible(path)
}
