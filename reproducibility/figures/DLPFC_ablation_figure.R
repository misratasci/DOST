# DLPFC12 ablation figures for DOST: gene-set size, embedding dimension,
# initialization and the spatial weight lambda.
#
# Reads the files the ablation scripts in reproducibility/DOST/ write, so run
# those first, then set dir.output below to the folder they wrote to.
#
# Inputs (all from dir.output). The per-slice files are indexed by slice
# position 1-12, not by slice ID:
#   DOST_DLPFC_ablation_nGenes.R  DOST_DLPFC_ablation_nGenes_HVG.RData  ari_means
#                                 DOST_DLPFC_ablation_nGenes_SVG.RData  ari_means
#   DOST_DLPFC_ablation_embdim.R  DOST_DLPFC<index>_embdim.RData        results, aris
#   DOST_DLPFC_ablation_init.R    DOST_DLPFC<index>_init_<init>.RData    results, ari
#   DOST_DLPFC_ablation_max_iter.R
#                                 DOST_DLPFC_ablation_max_iter_<index>.RData
#                                                                       results
#   DOST_DLPFC_ablation_lambda.R  DOST_DLPFC_ablation_lambda_<index>.RData
#                                                                       results
#   DOST_DLPFC_ablation_lambda_sensitivity.R
#                                 DOST_DLPFC<index>_lambda_noref.RData
#                                 DOST_DLPFC<index>_lambda_ref.RData     results, aris
#
# Outputs written to dir.figures:
#   nGenes.pdf                                 mean ARI against HVG and SVG count
#   DLPFC_ablation_embedding_dim_<index>.pdf   clustering and UMAP per embedding dim
#   DLPFC_initialization_<index>.pdf           clustering and UMAP per initialization
#   DLPFC_ablation_max_iter_<index>.pdf        clustering and UMAP per iteration budget
#   DLPFC_ablation_lambda_<index>.pdf          clustering and UMAP per lambda
#   DLPFC_lambda_sensitivity.pdf               ARI against lambda, faceted by slice

source("reproducibility/figures/load_data.R")
source("reproducibility/figures/figure_utils.R")

library(dplyr)

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to the folder the ablation scripts wrote their results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <- "path/to/figures"

# The embedding-dimension, initialization and max-iterations panels are each
# drawn for slices 1 and 9. Choose which one to use for each panel here.
embdim_slice_index <- 1
init_slice_index <- 1
max_iter_slice_index <- 1
lambda_panel_slice_index <- 1

# Slices the lambda sensitivity figure facets over
lambda_slice_indices <- 1:12
# ---------------------------------------------------------------------------

set.seed(1999)

slices <- c(151507:151510, 151669:151676)

layer_levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Gene-set size
# ---------------------------------------------------------------------------

nGenes_plot <- c("2000", "3000", "4000", "5000")

hvg_means <- load_object(file.path(dir.output, "DOST_DLPFC_ablation_nGenes_HVG.RData"),
                         "ari_means")
svg_means <- load_object(file.path(dir.output, "DOST_DLPFC_ablation_nGenes_SVG.RData"),
                         "ari_means")

df_hvg <- data.frame(nGenes = as.numeric(nGenes_plot),
                     ARI = as.numeric(hvg_means[nGenes_plot]),
                     Type = "HVG")
df_svg <- data.frame(nGenes = as.numeric(nGenes_plot),
                     ARI = as.numeric(svg_means[nGenes_plot]),
                     Type = "SVG")

plot_ablation <- function(data, title, x_label) {
  ggplot(data, aes(x = nGenes, y = ARI)) +
    geom_line() +
    geom_point(size = 3) +
    scale_x_continuous(breaks = data$nGenes) +
    labs(title = title, x = x_label, y = "Mean ARI") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p1 <- plot_ablation(df_hvg, "HVG", "Number of HVGs")
p2 <- plot_ablation(df_svg, "SVG", "Number of SVGs") + labs(y = NULL)

save_cropped(file.path(dir.figures, "nGenes.pdf"), p1 + p2,
             width = 8, height = 3.5)

# ---------------------------------------------------------------------------
# Embedding dimension
# ---------------------------------------------------------------------------

sample <- load_DLPFC_sample(slices[embdim_slice_index], dir.input)
sample$layers <- factor(sample$layers, levels = layer_levels)
gt <- sample@meta.data$layers

embdim_results <- load_object(
  file.path(dir.output, paste0("DOST_DLPFC", embdim_slice_index, "_embdim.RData")),
  "results")
embdims <- names(embdim_results)

spatial_plots <- list()
umap_plots <- list()
for (embedding_dim in embdims) {
  sample@meta.data$dost <- as.factor(embdim_results[[embedding_dim]]$labels)
  spatial_plots[[embedding_dim]] <-
    make_spatial_plot(sample, paste0("Embedding Dim. = ", embedding_dim),
                      "dost", "layers")
  umap_plots[[embedding_dim]] <-
    make_umap_plot(embdim_results[[embedding_dim]]$Z, gt, title = "",
                   trim_quant = 0.01)
}

embdim_grid <- wrap_plots(spatial_plots, ncol = length(embdims)) /
  wrap_plots(umap_plots, ncol = length(embdims))

save_cropped(file.path(dir.figures,
                       paste0("DLPFC_ablation_embedding_dim_",
                              embdim_slice_index, ".pdf")),
             embdim_grid,
             width = 2.5 * length(embdims), height = 6)

# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------

sample <- load_DLPFC_sample(slices[init_slice_index], dir.input)
sample$layers <- factor(sample$layers, levels = layer_levels)
gt <- sample@meta.data$layers

init_types <- c("random", "pca", "mds")
init_titles <- c(random = "Random", pca = "PCA", mds = "MDS")

init_columns <- list()
for (init in init_types) {
  init_file <- file.path(dir.output,
                         paste0("DOST_DLPFC", init_slice_index, "_init_", init,
                                ".RData"))
  if (!file.exists(init_file)) {
    message("skipping initialization '", init, "': ", basename(init_file), " not found")
    next
  }
  results <- load_object(init_file, "results")
  ari <- load_object(init_file, "ari")

  sample@meta.data$dost <- as.factor(results$labels)
  p_spatial <- Seurat::SpatialDimPlot(sample, group.by = "dost",
                                      pt.size.factor = 2.5, alpha = 1) +
    labs(caption = paste0("ARI = ", formatC(ari, digits = 3, format = "f"))) +
    Seurat::NoLegend() +
    theme(plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 10),
          plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0.5)) +
    ggtitle(init_titles[[init]])

  p_umap <- make_umap_plot(results$Z, gt, title = "", trim_quant = 0.01) +
    theme(plot.title = element_blank())

  init_columns[[init]] <- p_spatial / p_umap
}

init_grid <- wrap_plots(init_columns, ncol = length(init_columns))
save_cropped(file.path(dir.figures,
                       paste0("DLPFC_initialization_", init_slice_index, ".pdf")),
             init_grid,
             width = 3 * length(init_columns), height = 6)

# ---------------------------------------------------------------------------
# Max. iterations
# ---------------------------------------------------------------------------

max_iter_file <- file.path(dir.output,
                           paste0("DOST_DLPFC_ablation_max_iter_",
                                  max_iter_slice_index, ".RData"))

if (!file.exists(max_iter_file)) {
  message("skipping the max-iterations panel: ", basename(max_iter_file),
          " not found")
} else {
  sample <- load_DLPFC_sample(slices[max_iter_slice_index], dir.input)
  sample$layers <- factor(sample$layers, levels = layer_levels)
  gt <- sample@meta.data$layers

  max_iter_results <- load_object(max_iter_file, "results")
  max_iters <- names(max_iter_results)

  spatial_plots <- list()
  umap_plots <- list()
  for (max_iterations in max_iters) {
    sample@meta.data$dost <- as.factor(max_iter_results[[max_iterations]]$labels)
    # At 0 the optimizer returns its starting point, so that column is the raw MDS embedding
    title <- if (max_iterations == "0") {
      "0 iterations (MDS)"
    } else {
      paste0(max_iterations, " iterations")
    }
    spatial_plots[[max_iterations]] <-
      make_spatial_plot(sample, title, "dost", "layers")
    umap_plots[[max_iterations]] <-
      make_umap_plot(max_iter_results[[max_iterations]]$Z, gt, title = "",
                     trim_quant = 0.01)
  }

  max_iter_grid <- wrap_plots(spatial_plots, ncol = length(max_iters)) /
    wrap_plots(umap_plots, ncol = length(max_iters))

  save_cropped(file.path(dir.figures,
                         paste0("DLPFC_ablation_max_iter_",
                                max_iter_slice_index, ".pdf")),
               max_iter_grid,
               width = 2.5 * length(max_iters), height = 6)
}

# ---------------------------------------------------------------------------
# Lambda, fine grid around the default
# ---------------------------------------------------------------------------

lambda_panel_file <- file.path(dir.output,
                               paste0("DOST_DLPFC_ablation_lambda_",
                                      lambda_panel_slice_index, ".RData"))

if (!file.exists(lambda_panel_file)) {
  message("skipping the lambda panel: ", basename(lambda_panel_file), " not found")
} else {
  sample <- load_DLPFC_sample(slices[lambda_panel_slice_index], dir.input)
  sample$layers <- factor(sample$layers, levels = layer_levels)
  gt <- sample@meta.data$layers

  lambda_panel_results <- load_object(lambda_panel_file, "results")
  panel_lambdas <- names(lambda_panel_results)

  spatial_plots <- list()
  umap_plots <- list()
  for (lambda in panel_lambdas) {
    sample@meta.data$dost <- as.factor(lambda_panel_results[[lambda]]$labels)
    spatial_plots[[lambda]] <-
      make_spatial_plot(sample, paste0("Lambda = ", lambda), "dost", "layers")
    umap_plots[[lambda]] <-
      make_umap_plot(lambda_panel_results[[lambda]]$Z, gt, title = "",
                     trim_quant = 0.01)
  }

  lambda_panel_grid <- wrap_plots(spatial_plots, ncol = length(panel_lambdas)) /
    wrap_plots(umap_plots, ncol = length(panel_lambdas))

  save_cropped(file.path(dir.figures,
                         paste0("DLPFC_ablation_lambda_",
                                lambda_panel_slice_index, ".pdf")),
               lambda_panel_grid,
               width = 2.5 * length(panel_lambdas), height = 6)
}

# ---------------------------------------------------------------------------
# Lambda sensitivity
# ---------------------------------------------------------------------------

# One row per (slice, lambda) for one of the two refinement settings
load_lambda_results <- function(slice_index, is_refined) {
  suffix <- if (is_refined) "ref" else "noref"
  path <- file.path(dir.output,
                    paste0("DOST_DLPFC", slice_index, "_lambda_", suffix,
                           ".RData"))
  if (!file.exists(path)) return(NULL)
  aris <- load_object(path, "aris")
  data.frame(Lambda = as.numeric(names(aris)),
             ARI = as.numeric(aris),
             Slice_ID = as.character(slices[slice_index]),
             Method = if (is_refined) "DOST" else "DOST w/o Refinement",
             stringsAsFactors = FALSE)
}

lambda_results <- bind_rows(lapply(lambda_slice_indices, function(i) {
  rbind(load_lambda_results(i, FALSE), load_lambda_results(i, TRUE))
}))

if (nrow(lambda_results) == 0) {
  message("no lambda ablation files found in ", dir.output,
          ", skipping DLPFC_lambda_sensitivity.pdf")
} else {
  lambda_results$Slice_ID <- factor(lambda_results$Slice_ID,
                                    levels = as.character(slices))

  lambda_plot <- ggplot(lambda_results,
                        aes(x = Lambda, y = ARI, color = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    facet_wrap(~Slice_ID, ncol = 4, scales = "fixed") +
    labs(x = "Lambda (Spatial Regularization)",
         y = "Adjusted Rand Index (ARI)", color = NULL) +
    scale_color_manual(values = c("DOST w/o Refinement" = "#377EB8",
                                  "DOST" = "#E41A1C")) +
    theme_classic() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line.x = element_blank(),
          legend.position = "bottom")

  save_cropped(file.path(dir.figures, "DLPFC_lambda_sensitivity.pdf"),
               lambda_plot, width = 11, height = 9)
}
