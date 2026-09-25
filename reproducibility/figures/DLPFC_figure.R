# DLPFC12 figures: the slice grid and the ARI violin plot.
#
# Reads the files the scripts in reproducibility/<method>/ write, so run those
# first with their dir.output all pointing at the same folder, then set
# dir.output below to that folder. The Python methods run one slice per call,
# so all twelve slices have to be run before this script.
#
# Inputs (all from dir.output):
#   BASS       DLPFC_<slice>_BASS_labels.RData        zlabels
#   ADEPT      DLPFC_<slice>_adept.txt                one label per line
#   GraphST    graphst_slice_<slice>.csv              domain
#   stCluster  stcluster_slice_<slice>.csv            mclust
#   STAGATE    stagate_slice_<slice>.csv              mclust
#   SpaGCN     spagcn_slice_<slice>.csv               refined_pred
#   BANKSY     DLPFC_BANKSY_slice_<slice>_labels.RData  labels
#   DR.SC      DLPFC_DRSC_labels_<slice>.RData        labels
#   DOST       DOST_DLPFC.RData                       results
#              DOST_DLPFC_noref.RData                 results
#
# Outputs written to dir.figures:
#   DLPFC_slices<selected>.pdf   ground truth and every method, selected slices
#   DLPFC_violinplot.pdf         ARI over the 12 slices, one violin per method

source("reproducibility/figures/load_data.R")
source("reproducibility/figures/figure_utils.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/DLPFC12"

# Change this path to the folder the method scripts wrote their results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <- "path/to/figures"

# Slices shown as columns in the grid figure (indices into `slices`).
selected_indices <- c(1, 8, 9, 12)
# ---------------------------------------------------------------------------

set.seed(1999)

# Slice IDs
slices <- c(151507:151510, 151669:151676)
slice_indices <- 1:12

layer_levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")
layer_cols <- setNames(hue_pal()(length(layer_levels)), layer_levels)

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

# DOST saves all 12 slices in one file per refinement setting
dost_results <- load_object(file.path(dir.output, "DOST_DLPFC.RData"), "results")
dost_noref_results <- load_object(file.path(dir.output, "DOST_DLPFC_noref.RData"),
                                  "results")

# ---------------------------------------------------------------------------
# Data and ground truth
# ---------------------------------------------------------------------------

samples <- list()
for (i in slice_indices) {
  samples[[i]] <- load_DLPFC_sample(slices[i], dir.input)
  samples[[i]]$layers <- factor(samples[[i]]$layers, levels = layer_levels)
}

barcodes <- lapply(samples, colnames)

# ---------------------------------------------------------------------------
# Method labels, one metadata column per method
# ---------------------------------------------------------------------------

for (i in slice_indices) {
  ids <- barcodes[[i]]
  slice <- slices[i]

  zlabels <- load_object(file.path(dir.output,
                                   paste0("DLPFC_", slice, "_BASS_labels.RData")),
                         "zlabels")
  samples[[i]]@meta.data$bass <- zlabels[[1]]

  # ADEPT writes labels only, in the order its loader read the spots
  samples[[i]]@meta.data$adept <- as.factor(
    read_labels_txt(file.path(dir.output, paste0("DLPFC_", slice, "_adept.txt"))))

  samples[[i]]@meta.data$graphst <- as.factor(
    read_labels_csv(file.path(dir.output, paste0("graphst_slice_", slice, ".csv")),
                    "domain", ids))

  samples[[i]]@meta.data$stcluster <- as.factor(
    read_labels_csv(file.path(dir.output, paste0("stcluster_slice_", slice, ".csv")),
                    "mclust", ids))

  samples[[i]]@meta.data$stagate <- as.factor(
    read_labels_csv(file.path(dir.output, paste0("stagate_slice_", slice, ".csv")),
                    "mclust", ids))

  samples[[i]]@meta.data$spagcn <- as.factor(
    read_labels_csv(file.path(dir.output, paste0("spagcn_slice_", slice, ".csv")),
                    "refined_pred", ids))

  samples[[i]]@meta.data$banksy <- as.factor(
    load_object(file.path(dir.output,
                          paste0("DLPFC_BANKSY_slice_", slice, "_labels.RData")),
                "labels"))

  samples[[i]]@meta.data$drsc <- as.factor(
    load_object(file.path(dir.output,
                          paste0("DLPFC_DRSC_labels_", slice, ".RData")),
                "labels"))

  samples[[i]]@meta.data$dost <- as.factor(dost_results[[i]]$labels)

  samples[[i]]@meta.data$dost_noref <- as.factor(dost_noref_results[[i]]$labels)
}

# ---------------------------------------------------------------------------
# Slice grid
# ---------------------------------------------------------------------------

selected <- samples[selected_indices]

gt_col <- make_spatial_column(selected, "Ground Truth", "layers", "layers",
                              cols = layer_cols,
                              captions = paste0("Slice ", slices[selected_indices]))

method_columns <- list(
  c(title = "BASS",          column = "bass"),
  c(title = "ADEPT",         column = "adept"),
  c(title = "GraphST",       column = "graphst"),
  c(title = "stCluster",     column = "stcluster"),
  c(title = "STAGATE",       column = "stagate"),
  c(title = "SpaGCN",        column = "spagcn"),
  c(title = "BANKSY",        column = "banksy"),
  c(title = "DR.SC",         column = "drsc"),
  c(title = "DOST w/o ref.", column = "dost_noref"),
  c(title = "DOST",          column = "dost")
)

columns <- c(list(gt_col),
             lapply(method_columns, function(m) {
               make_spatial_column(selected, m[["title"]], m[["column"]], "layers")
             }))

grid <- wrap_plots(columns, ncol = length(columns))
save_cropped(file.path(dir.figures,
                       paste0("DLPFC_slices",
                              paste(selected_indices, collapse = "-"), ".pdf")),
             grid,
             width = 1.3 * length(columns),
             height = 1.5 * length(selected_indices))

# ---------------------------------------------------------------------------
# ARI violin plot over all 12 slices
# ---------------------------------------------------------------------------

aris <- sapply(method_columns, function(m) {
  sapply(1:length(slices), function(ind) {
    mclust::adjustedRandIndex(
      samples[[ind]]@meta.data$layers,
      samples[[ind]]@meta.data[[m[["column"]]]]
    )
  })
})

aris <- as.data.frame(aris, check.names = FALSE)
rownames(aris) <- slices
colnames(aris) <- sapply(method_columns, function(m) m[["title"]])

violin_plot <- make_ari_violin(aris)
save_cropped(file.path(dir.figures, "DLPFC_violinplot.pdf"),
             violin_plot, width = 8, height = 4)
