# mHypothalamus ablation figure: what the spatial weight lambda does.
#
# One page per section, two blocks. The top block sweeps lambda with R set to
# the number of cell classes and scores against the cell-class annotation; the
# bottom block sweeps lambda with R set to the number of domains and scores
# against the domain annotation. Together they show lambda moving DOST from
# cell typing to domain segmentation.
#
# Reads the files DOST_mHypothalamus_ablation_lambda.R writes, so run that
# first, then set dir.output below to the folder it wrote to.
#
# Inputs (all from dir.output):
#   mHypothalamus_DOST_ablation_cell_type_sheet_<sheet>.RData      results
#   mHypothalamus_DOST_ablation_spatial_domain_sheet_<sheet>.RData results
#
# Outputs written to dir.figures:
#   mHypothalamus_ablation_<sheet>.pdf   one per section

source("reproducibility/figures/load_data.R")
source("reproducibility/figures/figure_utils.R")

library(grid)

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mHypothalamus"

# Change this path to the folder the ablation script wrote its results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <- "path/to/figures"
# ---------------------------------------------------------------------------

set.seed(1999)

sheet <- '-0.24'

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

make_lambda_block <- function(results, coords, gtlabels, lambdas, sheet,
                              legend_ncol = 2, legend_spacer = FALSE, colors = c(0, 180)) {
  gt_plot <- make_points_plot(coords, gtlabels,
                              title = "Ground Truth",
                              caption = paste0("Bregma ", sheet),
                              size = 1, legend = FALSE, colors = colors)
  gt_legend <- make_points_plot(coords, gtlabels,
                                title = "Ground Truth",
                                caption = paste0("Bregma ", sheet),
                                size = 1, legend = TRUE,
                                legend_ncol = legend_ncol, colors = colors)
  gt_legend <- wrap_elements(cowplot::get_legend(gt_legend))

  spatial_plots <- list()
  umap_plots <- list()
  for (lambda in lambdas) {
    key <- as.character(lambda)
    labels <- as.factor(results[[key]]$labels)
    ari <- mclust::adjustedRandIndex(labels, gtlabels)
    spatial_plots[[key]] <- make_points_plot(
      coords, labels,
      title = paste0("Lambda = ", lambda),
      caption = paste0("ARI=", formatC(ari, digits = 2, format = "f")),
      size = 1, colors = colors)
    umap_plots[[key]] <- make_umap_plot(results[[key]]$Z, gtlabels, title = "",
                                        trim_quant = 0.01, colors = colors)
  }
  if (legend_spacer) gt_legend <- gt_legend / plot_spacer()
  spatial_row <- (gt_plot | wrap_plots(spatial_plots, ncol = length(lambdas))) +
    plot_layout(widths = c(0.2, 1))
  umap_row <- (gt_legend | wrap_plots(umap_plots, ncol = length(lambdas))) +
    plot_layout(widths = c(0.2, 1))
  spatial_row / umap_row
}

cell_file <- file.path(dir.output,
                       paste0("mHypothalamus_DOST_ablation_cell_type_sheet_",
                              sheet, ".RData"))
domain_file <- file.path(dir.output,
                         paste0("mHypothalamus_DOST_ablation_spatial_domain_sheet_",
                                sheet, ".RData"))

celltype_res <- load_object(cell_file, "results")
spatialdomain_res <- load_object(domain_file, "results")

lambdas <- as.numeric(names(spatialdomain_res))

info <- load_mHypothalamus_info(dir.input, sheet)

cell_block <- make_lambda_block(celltype_res, info$coords,
                                info$cell_class, lambdas, sheet, colors = c(0, 180))
domain_block <- make_lambda_block(spatialdomain_res, info$coords,
                                  info$domains, lambdas, sheet,
                                  legend_ncol = 3, legend_spacer = TRUE, colors = c(180, 360))

title_cell <- wrap_elements(full = textGrob(
  "Cell types", rot = 90, gp = gpar(fontsize = 14, fontface = "bold")))
title_domain <- wrap_elements(full = textGrob(
  "Spatial domains", rot = 90, gp = gpar(fontsize = 14, fontface = "bold")))

row1 <- (title_cell | cell_block) + plot_layout(widths = c(0.02, 0.98))
row2 <- (title_domain | domain_block) + plot_layout(widths = c(0.02, 0.98))

save_cropped(file.path(dir.figures,
                       paste0("mHypothalamus_ablation_", sheet, ".pdf")),
             row1 / row2,
             width = 2.4 * (length(lambdas) + 1), height = 12)
