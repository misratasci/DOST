source("R/funcs_.R")

set.seed(1999)
section <- "MA"

sample <- load_mMAMP_sample("data/mMAMP", section)
R <- length(unique(sample@meta.data$ground_truth))
X <- Seurat::GetAssayData(sample, layer = "counts")
coords <- Seurat::GetTissueCoordinates(sample, scale = 'hires')[1:2]
results <- DOST(X, coords, R, selected_genes = "HVG", nGenes = 3000,
                neighborhood_threshold = 1, embedding_dim = 20,
                lambda = 0.03, lr = 16, max_iterations = 20,
                loss_tol = 1e-5, eps = 1e-20,
                refinement = FALSE)

plot_sample(sample, results$labels, sample@meta.data$ground_truth)

plot_umap(results$Z, sample@meta.data$ground_truth)

save(results, file = "RData/results_mMAMP_DOST.RData")
