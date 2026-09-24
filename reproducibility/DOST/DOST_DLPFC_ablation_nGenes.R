# Ablation: how many genes DOST needs, on DLPFC12.
#
# Sweeps the gene-set size over all 12 slices, once with highly variable genes
# and once with spatially variable genes, and records the mean ARI at each size.
#
# SVG selection needs the SPARK package:
#   devtools::install_github('xzhoulab/SPARK')
#
# Outputs written to dir.output:
#   DOST_DLPFC_ablation_nGenes_HVG.RData   aris, ari_means
#   DOST_DLPFC_ablation_nGenes_SVG.RData   aris, ari_means
#
# aris is a slices x gene-counts matrix, ari_means its column means, named by
# gene count.

library(DOST)
source("reproducibility/DOST/load_data.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir.output <- "path/to/output"

# Gene-set sizes to sweep and which selection methods to run
nGenes_grid <- c(2000, 3000, 4000, 5000)
gene_sets <- c("HVG", "SVG")

# Slice indices to average over (1 to 12)
slice_indices <- 1:12
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

for (gene_set in gene_sets) {

  aris <- matrix(NA_real_,
                 nrow = length(slice_indices), ncol = length(nGenes_grid),
                 dimnames = list(as.character(slices[slice_indices]),
                                 as.character(nGenes_grid)))

  for (slice_index in slice_indices) {

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

    for (nGenes in nGenes_grid) {
      results <- DOST(X, coords, R,
                      selected_genes = gene_set,
                      nGenes = nGenes)

      ari <- mclust::adjustedRandIndex(results$labels, sample@meta.data$layers)
      aris[as.character(slices[slice_index]), as.character(nGenes)] <- ari

      cat(gene_set, " slice:", slices[slice_index], " nGenes:", nGenes,
          " ARI:", ari, "\n")
    }
  }

  ari_means <- colMeans(aris, na.rm = TRUE)
  print(round(ari_means, 3))

  save(aris, ari_means,
       file = file.path(dir.output,
                        paste0("DOST_DLPFC_ablation_nGenes_", gene_set, ".RData")))
}
