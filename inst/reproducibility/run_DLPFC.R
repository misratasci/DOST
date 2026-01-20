source("R/funcs_.R")
library(umap)
library(scales)
library(Seurat)

dir.input <- "~/Desktop/ST-research/data/DLPFC12"
slices <- c(151507:151510, 151669:151676)
eps <- 1e-20

nGenes <- 1000
neighborhood_threshold <- 1
embedding_dim <- 20
selected_genes <- "SVG" # "HVG" or "SVG"
lambda <- 0.03
lr <- 16
max_iterations <- 20
loss_tol <- 1e-5
refinement <- TRUE
refine_k <- 30

aris <- c()
results <- list()
for (slice_index in 1:12) {
  cat(paste0("\nSlice ", slice_index, "\n"))
  if (slice_index > 4 & slice_index < 9) {
    R <- 5
  } else { 
    R <- 7
  }
  sample <- load_DLPFC_sample(slices[slice_index], dir.input)
  X <- GetAssayData(sample, layer = "counts")
  coords <- GetTissueCoordinates(sample, scale = 'hires')[1:2]
  results[[slice_index]] <- DOST(X, coords, R, selected_genes = selected_genes,
                                 nGenes = nGenes, neighborhood_threshold = neighborhood_threshold,
                                 embedding_dim = embedding_dim, lambda = lambda,
                                 lr = lr, max_iterations = max_iterations,
                                 loss_tol = loss_tol, eps = eps,
                                 refinement = refinement, refine_k = refine_k)
  aris <- c(aris, mclust::adjustedRandIndex(results[[slice_index]]$labels, sample@meta.data$layers))
}
save(results, aris, file = paste0("RData/DOST_DLPFC_nGenes1000_SVG.RData"))

results <- list()
for (max_iterations in seq(0, 60, 20)) {
  slice_index <- 9
  cat(paste0("\nSlice ", slice_index, "\n"))
  if (slice_index > 4 & slice_index < 9) {
    R <- 5
  } else { 
    R <- 7
  }
  sample <- load_DLPFC_sample(slices[slice_index], dir.input)
  X <- GetAssayData(sample, layer = "counts")
  coords <- GetTissueCoordinates(sample, scale = 'hires')[1:2]
  results[[as.character(max_iterations)]] <- DOST(X, coords, R, selected_genes = "HVG",
                  nGenes = 3000, neighborhood_threshold = 1, embedding_dim = 20,
                  lambda = 0.03, lr = 16, max_iterations = max_iterations,
                  loss_tol = 1e-5, eps = 1e-20,
                  refinement = TRUE, refine_k = 30)
}
save(results, file = paste0("RData/DOST_DLPFC_ablation_max_iter_9_20.RData"))

results <- list()
for (lambda in seq(0, 0.09, 0.03)) {
  slice_index <- 1
  cat(paste0("\nSlice ", slice_index, "\n"))
  if (slice_index > 4 & slice_index < 9) {
    R <- 5
  } else { 
    R <- 7
  }
  sample <- load_DLPFC_sample(slices[slice_index], dir.input)
  X <- GetAssayData(sample, layer = "counts")
  coords <- GetTissueCoordinates(sample, scale = 'hires')[1:2]
  results[[as.character(lambda)]] <- DOST(X, coords, R, selected_genes = "HVG",
                                                  nGenes = 3000, neighborhood_threshold = 1, embedding_dim = embedding_dim,
                                                  lambda = lambda, lr = 16, max_iterations = max_iterations,
                                                  loss_tol = 1e-5, eps = 1e-20,
                                                  refinement = TRUE, refine_k = 30)
}
save(results, file = paste0("RData/DOST_DLPFC_ablation_lambda_1.RData"))

results <- list()
aris <- matrix(0, nrow = 12, ncol = 5)
ari_means <- c()
i <- 0
for (nGenes in c(1000, 2000, 3000, 4000, 5000)) {
  for (slice_index in 1:12) {
    cat(paste0("\nSlice ", slice_index, "\n"))
    if (slice_index > 4 & slice_index < 9) {
      R <- 5
    } else { 
      R <- 7
    }
    sample <- load_DLPFC_sample(slices[slice_index], dir.input)
    X <- GetAssayData(sample, layer = "counts")
    coords <- GetTissueCoordinates(sample, scale = 'hires')[1:2]
    results[[slice_index]] <- DOST(X, coords, R, selected_genes = "SVG",
                                  nGenes = nGenes, neighborhood_threshold = 1, embedding_dim = 20,
                                  lambda = 0.03, lr = 16, max_iterations = 20,
                                  loss_tol = 1e-5, eps = 1e-20,
                                  refinement = TRUE, refine_k = 30)
    aris[slice_index, i] <- mclust::adjustedRandIndex(results[[slice_index]]$labels, sample@meta.data$layers)
  }
  ari_means <- c(ari_means, mean(aris[, i]))
  i <- i + 1
}
save(aris, ari_means, file = paste0("RData/DOST_DLPFC_ablation_nGenes_SVG.RData"))

library(tidyverse)
library(magick)
library(patchwork)

plots <- list()
for (i in 1:20) {
  umap_plot <- plot_umap(scale(results$Zs[[i]]), sample@meta.data$layers)
  plot <- plot_sample(sample,
                      cluster_embedding(scale(results$Zs[[i]]), 7)$classification,
                      sample@meta.data$layers)
  plots[[i]] <- wrap_plots(umap_plot, plot, ncol = 2)
}
png_files <- c()
for (i in 1:20) {
  ggplot2::ggsave(plots[[i]], file = paste0("iter_", i, ".png"), width = 8, height = 4)
  png_files <- c(png_files, paste0("iter_", i, ".png"))
}

### create a GIF file from all the plots
png_files %>%
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps = 1) %>% # animates
  image_write("All_plots.gif")

for (i in 1:12) {
  sample <- load_DLPFC_sample(slices[i], dir.input)
  write.csv(as.numeric(factor(sample@meta.data$layers)), file = paste0("DOST_slice", i, "_gt.csv"), row.names = FALSE)
}
