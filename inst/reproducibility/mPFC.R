
library(DOST)
set.seed(1999)

# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download mHypothalamus data from https://zenodo.org/records/10698914
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/"

# Load data
filename = paste0(dir.input, 'starmap_mpfc_starmap_cnts.xlsx')
infoname = paste0(dir.input, 'starmap_mpfc_starmap_info.xlsx')

# Read sheet IDs
sheets <- readxl::excel_sheets(filename)

# Run DOST for a selected sheet
sheet <- sheets[1]

# Load data
cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
row.names(cnts) <- cnts[,"...1"]
cnts <- cnts[ -c(1) ]
xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
row.names(xys) <- xys[,"...1"]
gtlabels <- xys$z
xys <- xys[-c(1)]
xys <- xys[-c(-2:-1)]
xys <- xys[,c(2,1)]
R <- length(unique(gtlabels))

# Run DOST
results <- DOST(cnts, xys, R,
                neighborhood_threshold = 1,
                embedding_dim = 10,
                lambda = 0.15,
                refinement = TRUE)

# Plot results
plot_DOST_spatial(xys, results$labels)
plot_DOST_umap(results$Z, gtlabels)

# Calculate ARI
mclust::adjustedRandIndex(results$labels, gtlabels)
