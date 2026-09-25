# mHypothalamus (MERFISH mouse hypothalamus) figures: the section grid and the
# ARI plots over the five annotated sections.
#
# Reads the files the scripts in reproducibility/<method>/ write, so run those
# first with their dir.output all pointing at the same folder, then set
# dir.output below to that folder. The Python methods run one section per call,
# so all five sections have to be run before this script.
#
# Inputs (all from dir.output):
#   BASS       BASS_mHypothalamus.RData                       results
#   ADEPT      mHypothalamus_<sheet>_adept.txt                one label per line
#   GraphST    graphst_mHypothalamus_<sheet>.csv              domain
#   stCluster  stcluster_mhypo_<sheet>.csv                    mclust
#   STAGATE    stagate_mHypothalamus_<sheet>.csv              mclust
#   SpaGCN     spagcn_mHypothalamus_<sheet>.csv               refined_pred
#   BANKSY     mHypothalamus_BANKSY_sheet_<sheet>_labels.RData  labels
#   DR.SC      mHypothalamus_DRSC_sheet_<sheet>_labels.RData    labels
#   DOST       mHypothalamus_DOST_sheet_<sheet>_labels.RData     labels
#
# Outputs written to dir.figures:
#   mHypo_slices.pdf         ground truth and every method, five sections
#   mHypo_violinplot.pdf     ARI over the five sections, one violin per method

source("reproducibility/figures/load_data.R")
source("reproducibility/figures/figure_utils.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mHypothalamus"

# Change this path to the folder the method scripts wrote their results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <- "path/to/figures"
# ---------------------------------------------------------------------------

set.seed(1999)

# Sheet IDs with domain annotations
sheets <- c('-0.04', '-0.09', '-0.14', '-0.19', '-0.24')

# Dot size in the section grid
dotsize <- 0.1

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Coordinates and ground truth
# ---------------------------------------------------------------------------

info <- lapply(sheets, function(sheet) load_mHypothalamus_info(dir.input, sheet))
names(info) <- sheets

coords <- lapply(info, function(x) x$coords)
cell_ids <- lapply(info, function(x) x$cell_ids)

gtlabels <- lapply(info, function(x) {
  l <- factor(x$domains)
})

# ---------------------------------------------------------------------------
# Method labels, one list of five label vectors per method
# ---------------------------------------------------------------------------

bass_results <- load_object(file.path(dir.output, "BASS_mHypothalamus.RData"), "results")
bass_labels <- lapply(sheets, function(sheet) as.factor(bass_results[[sheet]]$z[[1]]))

adept_labels <- lapply(sheets, function(sheet) {
  as.factor(read_labels_txt(file.path(dir.output,
                                      paste0("mHypothalamus_", sheet, "_adept.txt"))))
})

graphst_labels <- lapply(sheets, function(sheet) {
  as.factor(read_labels_csv(file.path(dir.output,
                                      paste0("graphst_mHypothalamus_", sheet, ".csv")),
                            "domain", cell_ids[[sheet]]))
})

stcluster_labels <- lapply(sheets, function(sheet) {
  as.factor(read_labels_csv(file.path(dir.output,
                                      paste0("stcluster_mhypo_", sheet, ".csv")),
                            "mclust", cell_ids[[sheet]]))
})

stagate_labels <- lapply(sheets, function(sheet) {
  as.factor(read_labels_csv(file.path(dir.output,
                                      paste0("stagate_mHypothalamus_", sheet, ".csv")),
                            "mclust", cell_ids[[sheet]]))
})

spagcn_labels <- lapply(sheets, function(sheet) {
  as.factor(read_labels_csv(file.path(dir.output,
                                      paste0("spagcn_mHypothalamus_", sheet, ".csv")),
                            "refined_pred", cell_ids[[sheet]]))
})

banksy_labels <- lapply(sheets, function(sheet) {
  as.factor(load_object(file.path(dir.output,
                                  paste0("mHypothalamus_BANKSY_sheet_", sheet,
                                         "_labels.RData")), "labels"))
})

drsc_labels <- lapply(sheets, function(sheet) {
  as.factor(load_object(file.path(dir.output,
                                  paste0("mHypothalamus_DRSC_sheet_", sheet,
                                         "_labels.RData")), "labels"))
})

dost_labels <- lapply(sheets, function(sheet) {
  as.factor(load_object(file.path(dir.output,
                                  paste0("mHypothalamus_DOST_sheet_", sheet,
                                         "_labels.RData")), "labels"))
})

# ---------------------------------------------------------------------------
# Section grid
# ---------------------------------------------------------------------------

#sheets <- c('-0.04','-0.19') for figure in main manuscript with two selected sheets
annot_plots <- make_points_column(coords, gtlabels, gtlabels, sheets,
                                  "Ground Truth", size = dotsize,
                                  captions = paste0("Bregma ", sheets))
manual_annot_col <- wrap_plots(annot_plots, ncol = 1, guides = "collect")

method_labels <- list(
  "BASS"      = bass_labels,
  "ADEPT"     = adept_labels,
  "GraphST"   = graphst_labels,
  "stCluster" = stcluster_labels,
  "STAGATE"   = stagate_labels,
  "SpaGCN"    = spagcn_labels,
  "BANKSY"    = banksy_labels,
  "DR.SC"     = drsc_labels,
  "DOST"      = dost_labels
)

method_cols <- lapply(names(method_labels), function(method) {
  wrap_plots(make_points_column(coords, method_labels[[method]], gtlabels, sheets,
                                method, size = dotsize),
             ncol = 1)
})

columns <- c(list(manual_annot_col), method_cols)
grid <- wrap_plots(columns, ncol = length(columns))
save_cropped(file.path(dir.figures, "mHypo_slices.pdf"), grid,
             width = 1.3 * length(columns), height = 7.5)
# height = 3.0 for figure in main manuscript with 2 sheets

# ---------------------------------------------------------------------------
# ARI plots
# ---------------------------------------------------------------------------

aris <- sapply(method_labels, function(labels_by_sheet) {
  sapply(1:length(sheets), function(ind) {
    mclust::adjustedRandIndex(
      gtlabels[[ind]],
      labels_by_sheet[[ind]]
    )
  })
})

aris <- as.data.frame(aris, check.names = FALSE)
rownames(aris) <- sheets

violin_plot <- make_ari_violin(aris, show_mean = TRUE)
save_cropped(file.path(dir.figures, "mHypo_violinplot.pdf"), violin_plot,
             width = 6.5, height = 2.5)
