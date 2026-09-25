# BC (human breast cancer, Visium section1) figure.
#
# Reads the files the scripts in reproducibility/<method>/ write, so run those
# first with their dir.output all pointing at the same folder, then set
# dir.output below to that folder.
#
# Inputs (all from dir.output):
#   BASS       BC_BASS_labels.RData             zlabels
#   ADEPT      BC_adept.txt                     one label per line
#              adept_embeddings_BC.csv
#   GraphST    graphst_bc.csv                   domain
#              graphst_embeddings_bc.csv
#   stCluster  stcluster_bc_section1.csv        mclust
#              stcluster_embeddings_bc_section1.csv
#   STAGATE    stagate_bc.csv                   mclust
#              stagate_embeddings_bc.csv
#   SpaGCN     spagcn_bc.csv                    refined_pred
#              spagcn_embeddings_bc.csv
#   BANKSY     BC_BANKSY_labels.RData           labels
#              BC_BANKSY_embeddings.RData       emb
#   DR.SC      BC_DRSC_labels.RData             labels
#              BC_DRSC_embeddings.RData         emb
#   DOST       BC_DOST_labels.RData             labels
#              BC_DOST_embeddings.RData         emb
#
# Outputs written to dir.figures:
#   BC.pdf

source("reproducibility/figures/load_data.R")
source("reproducibility/figures/figure_utils.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/BC"

# Change this path to the folder the method scripts wrote their results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <- "path/to/figures"
# ---------------------------------------------------------------------------

set.seed(1999)

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

sample <- load_BC_sample(dir.input)
barcodes <- colnames(sample)
gt <- sample@meta.data$fine_annot_type

# ---------------------------------------------------------------------------
# Labels and embeddings
# ---------------------------------------------------------------------------

zlabels <- load_object(file.path(dir.output, "BC_BASS_labels.RData"), "zlabels")
sample@meta.data$bass <- zlabels[[1]]

sample@meta.data$adept <- as.factor(
  read_labels_txt(file.path(dir.output, "BC_adept.txt")))
adept_emb <- read_embedding_csv(file.path(dir.output, "adept_embeddings_BC.csv"))

graphst_csv <- read.csv(file.path(dir.output, "graphst_bc.csv"), row.names = 1)
sample@meta.data$graphst <- as.factor(
  align_to_ids(graphst_csv$domain, rownames(graphst_csv), barcodes, "graphst_bc.csv"))
graphst_emb <- read_embedding_csv(file.path(dir.output, "graphst_embeddings_bc.csv"),
                                  barcodes, rownames(graphst_csv))

stcluster_csv <- read.csv(file.path(dir.output, "stcluster_bc_section1.csv"), row.names = 1)
sample@meta.data$stcluster <- as.factor(
  align_to_ids(stcluster_csv$mclust, rownames(stcluster_csv), barcodes,
               "stcluster_bc_section1.csv"))
stcluster_emb <- read_embedding_csv(file.path(dir.output, "stcluster_embeddings_bc_section1.csv"),
                                    barcodes, rownames(stcluster_csv))

stagate_csv <- read.csv(file.path(dir.output, "stagate_bc.csv"), row.names = 1)
sample@meta.data$stagate <- as.factor(
  align_to_ids(stagate_csv$mclust, rownames(stagate_csv), barcodes, "stagate_bc.csv"))
stagate_emb <- read_embedding_csv(file.path(dir.output, "stagate_embeddings_bc.csv"),
                                  barcodes, rownames(stagate_csv))

spagcn_csv <- read.csv(file.path(dir.output, "spagcn_bc.csv"), row.names = 1)
sample@meta.data$spagcn <- as.factor(
  align_to_ids(spagcn_csv$refined_pred, rownames(spagcn_csv), barcodes, "spagcn_bc.csv"))
spagcn_emb <- read_embedding_csv(file.path(dir.output, "spagcn_embeddings_bc.csv"),
                                 barcodes, rownames(spagcn_csv))

sample@meta.data$banksy <- as.factor(
  load_object(file.path(dir.output, "BC_BANKSY_labels.RData"), "labels"))
banksy_emb <- load_object(file.path(dir.output, "BC_BANKSY_embeddings.RData"), "emb")

sample@meta.data$drsc <- as.factor(
  load_object(file.path(dir.output, "BC_DRSC_labels.RData"), "labels"))
drsc_emb <- load_object(file.path(dir.output, "BC_DRSC_embeddings.RData"), "emb")

sample@meta.data$dost <- as.factor(
  load_object(file.path(dir.output, "BC_DOST_labels.RData"), "labels"))
dost_emb <- load_object(file.path(dir.output, "BC_DOST_embeddings.RData"), "emb")

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

spatial_row <- list(
  make_spatial_plot(sample, "Ground Truth", "fine_annot_type", "fine_annot_type",
                    show_ari = FALSE),
  make_spatial_plot(sample, "BASS",      "bass",      "fine_annot_type"),
  make_spatial_plot(sample, "ADEPT",     "adept",     "fine_annot_type"),
  make_spatial_plot(sample, "GraphST",   "graphst",   "fine_annot_type"),
  make_spatial_plot(sample, "stCluster", "stcluster", "fine_annot_type"),
  make_spatial_plot(sample, "STAGATE",   "stagate",   "fine_annot_type"),
  make_spatial_plot(sample, "SpaGCN",    "spagcn",    "fine_annot_type"),
  make_spatial_plot(sample, "BANKSY",    "banksy",    "fine_annot_type"),
  make_spatial_plot(sample, "DR.SC",     "drsc",      "fine_annot_type"),
  make_spatial_plot(sample, "DOST",      "dost",      "fine_annot_type")
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

save_cropped(file.path(dir.figures, "BC.pdf"), final_plot,
             width = 1.35 * n_col / 2, height = 7.6)
