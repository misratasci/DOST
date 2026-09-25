# mMAMP (mouse brain anterior, Visium section MA) figure.
#
# Reads the files the scripts in reproducibility/<method>/ write, so run those
# first with their dir.output all pointing at the same folder, then set
# dir.output below to that folder.
#
# Inputs (all from dir.output):
#   BASS       mMAMP_BASS_labels.RData             zlabels
#   ADEPT      mMAMP_adept.txt                     one label per line
#              adept_embeddings_mMAMP.csv
#   GraphST    graphst_mMAMP.csv                   domain
#              graphst_embeddings_mMAMP.csv
#   stCluster  stcluster_mmamp_MA.csv              mclust
#              stcluster_embeddings_mmamp_MA.csv
#   STAGATE    stagate_mMAMP.csv                   mclust
#              stagate_embeddings_mMAMP.csv
#   SpaGCN     spagcn_mMAMP.csv                    refined_pred
#              spagcn_embeddings_mMAMP.csv
#   BANKSY     mMAMP_BANKSY_labels.RData           labels
#              mMAMP_BANKSY_embeddings.RData       emb
#   DR.SC      mMAMP_DRSC_labels.RData             labels
#              mMAMP_DRSC_embeddings.RData         emb
#   DOST       mMAMP_DOST_labels.RData             labels
#              mMAMP_DOST_embeddings.RData         emb
#
# Outputs written to dir.figures:
#   mMAMP.pdf

source("reproducibility/figures/load_data.R")
source("reproducibility/figures/figure_utils.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mMAMP"

# Change this path to the folder the method scripts wrote their results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <-  "path/to/figures"
# ---------------------------------------------------------------------------

set.seed(1999)

section <- "MA"

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

sample <- load_mMAMP_sample(dir.input, section)
barcodes <- colnames(sample)
gt <- sample@meta.data$ground_truth

# ---------------------------------------------------------------------------
# Labels and embeddings
# ---------------------------------------------------------------------------

# BASS: zlabels is a list with one entry per section, here always one
zlabels <- load_object(file.path(dir.output, "mMAMP_BASS_labels.RData"), "zlabels")
sample@meta.data$bass <- zlabels

# ADEPT writes labels only, in the order its loader read the spots
sample@meta.data$adept <- as.factor(read_labels_txt(file.path(dir.output, "mMAMP_adept.txt")))
adept_emb <- read_embedding_csv(file.path(dir.output, "adept_embeddings_mMAMP.csv"))

graphst_csv <- read.csv(file.path(dir.output, "graphst_mMAMP.csv"), row.names = 1)
sample@meta.data$graphst <- as.factor(
  align_to_ids(graphst_csv$domain, rownames(graphst_csv), barcodes, "graphst_mMAMP.csv"))
graphst_emb <- read_embedding_csv(file.path(dir.output, "graphst_embeddings_mMAMP.csv"),
                                  barcodes, rownames(graphst_csv))

stcluster_csv <- read.csv(file.path(dir.output, "stcluster_mmamp_MA.csv"), row.names = 1)
sample@meta.data$stcluster <- as.factor(
  align_to_ids(stcluster_csv$mclust, rownames(stcluster_csv), barcodes,
               "stcluster_mmamp_MA.csv"))
stcluster_emb <- read_embedding_csv(file.path(dir.output, "stcluster_embeddings_mmamp_MA.csv"),
                                    barcodes, rownames(stcluster_csv))

stagate_csv <- read.csv(file.path(dir.output, "stagate_mMAMP.csv"), row.names = 1)
sample@meta.data$stagate <- as.factor(
  align_to_ids(stagate_csv$mclust, rownames(stagate_csv), barcodes, "stagate_mMAMP.csv"))
stagate_emb <- read_embedding_csv(file.path(dir.output, "stagate_embeddings_mMAMP.csv"),
                                  barcodes, rownames(stagate_csv))

spagcn_csv <- read.csv(file.path(dir.output, "spagcn_mMAMP.csv"), row.names = 1)
sample@meta.data$spagcn <- as.factor(
  align_to_ids(spagcn_csv$refined_pred, rownames(spagcn_csv), barcodes, "spagcn_mMAMP.csv"))
spagcn_emb <- read_embedding_csv(file.path(dir.output, "spagcn_embeddings_mMAMP.csv"),
                                 barcodes, rownames(spagcn_csv))

sample@meta.data$banksy <- as.factor(
  load_object(file.path(dir.output, "mMAMP_BANKSY_labels.RData"), "labels"))
banksy_emb <- load_object(file.path(dir.output, "mMAMP_BANKSY_embeddings.RData"), "emb")

sample@meta.data$drsc <- as.factor(
  load_object(file.path(dir.output, "mMAMP_DRSC_labels.RData"), "labels"))
drsc_emb <- load_object(file.path(dir.output, "mMAMP_DRSC_embeddings.RData"), "emb")

sample@meta.data$dost <- as.factor(
  load_object(file.path(dir.output, "mMAMP_DOST_labels.RData"), "labels"))
dost_emb <- load_object(file.path(dir.output, "mMAMP_DOST_embeddings.RData"), "emb")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

spatial_row <- list(
  make_spatial_plot(sample, "Ground Truth", "ground_truth", "ground_truth",
                    show_ari = FALSE),
  make_spatial_plot(sample, "BASS",      "bass",      "ground_truth"),
  make_spatial_plot(sample, "ADEPT",     "adept",     "ground_truth"),
  make_spatial_plot(sample, "GraphST",   "graphst",   "ground_truth"),
  make_spatial_plot(sample, "stCluster", "stcluster", "ground_truth"),
  make_spatial_plot(sample, "STAGATE",   "stagate",   "ground_truth"),
  make_spatial_plot(sample, "SpaGCN",    "spagcn",    "ground_truth"),
  make_spatial_plot(sample, "BANKSY",    "banksy",    "ground_truth"),
  make_spatial_plot(sample, "DR.SC",     "drsc",      "ground_truth"),
  make_spatial_plot(sample, "DOST",      "dost",      "ground_truth")
)

umap_row <- list(
  make_umap_plot(adept_emb,     gt, "ADEPT"),
  make_umap_plot(graphst_emb,   gt, "GraphST"),
  make_umap_plot(stcluster_emb, gt, "stCluster"),
  make_umap_plot(stagate_emb,   gt, "STAGATE"),
  make_umap_plot(spagcn_emb,    gt, "SpaGCN"),
  make_umap_plot(banksy_emb,    gt, "BANKSY"),
  make_umap_plot(drsc_emb,      gt, "DR.SC"),
  make_umap_plot(dost_emb,      gt, "DOST")
)

n_col <- length(spatial_row)
final_plot <- wrap_plots(wrap_plots(spatial_row, ncol = n_col / 2),
                         wrap_plots(umap_row, ncol = (n_col - 2) / 2),
                         ncol = 1, heights = c(1, 1.2))

save_cropped(file.path(dir.figures, "mMAMP.pdf"), final_plot,
             width = 1.35 * n_col / 2, height = 7.6)
