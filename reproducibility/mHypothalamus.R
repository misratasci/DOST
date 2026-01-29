
library(DOST)
set.seed(1999)

# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download mHypothalamus data from https://zenodo.org/records/10698909
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/"

# Load data
filename = paste0(dir.input, 'MERFISH_Animal1_cnts.xlsx')
infoname = paste0(dir.input, 'MERFISH_Animal1_info.xlsx')

# Sheet IDs with labels
sheets <- c('-0.04', '-0.09', '-0.14', '-0.19', '-0.24')

# Run DOST for a selected sheet
sheet <- '-0.24'

# Load data
cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
row.names(cnts) <- cnts[,"...1"]
cnts <- cnts[-c(1)]
xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
row.names(xys) <- xys[,"...1"]
gtlabels <- xys$z
celllabels <- xys$Cell_class
xys <- xys[-c(1)]
xys <- xys[-c(-2:-1)]
xys <- xys[,c(2,1)]
R <- length(unique(gtlabels))

# Run DOST
results <- DOST(cnts, xys, R,
                neighborhood_threshold = 3,
                embedding_dim = 10,
                lambda = 0.2,
                refinement = FALSE)

# Plot results
plot_DOST_spatial(xys, results$labels)
plot_DOST_umap(results$Z, gtlabels)

# Calculate ARI
mclust::adjustedRandIndex(results$labels, gtlabels)

# Run DOST with 15 domains (number of unique cell type labels) and lambda = 0
results_cell <- DOST(cnts, xys, R = 15,
                     neighborhood_threshold = 3,
                     embedding_dim = 10,
                     lambda = 0,
                     refinement = FALSE)

# Plot results for cell types
plot_DOST_spatial(xys, results_cell$labels)
plot_DOST_umap(results_cell$Z, celllabels)

# Calculate ARI for cell types
mclust::adjustedRandIndex(results_cell$labels, celllabels)
