# Ablation: the spatial weight lambda, on mHypothalamus.
#
# Two sweeps per section over the same lambda grid:
#
#   cell_type       R = number of annotated cell classes, scored against the
#                   cell-class annotation
#   spatial_domain  R = number of annotated domains, scored against the domain
#                   annotation
#
# At lambda = 0 the spatial term is off and DOST behaves like a cell-type
# clusterer; as lambda grows the clusters turn into contiguous domains. Running
# both sweeps is what lets the figure show that in one page per section.
#
# We used the data links provided by Benchmark ST study:
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
# Download mHypothalamus data from https://zenodo.org/records/10698909
#
# Outputs written to dir.output:
#   mHypothalamus_DOST_ablation_cell_type_sheet_<sheet>.RData       results, aris
#   mHypothalamus_DOST_ablation_spatial_domain_sheet_<sheet>.RData  results, aris
#
# results is a list keyed by lambda as a string, each entry the full DOST
# return value (Z, losses, labels). aris is named the same way.

library(DOST)

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir.output <- "path/to/output/"

# Lambda grid to sweep
lambdas <- seq(0, 0.2, by = 0.05)

# Sections to run
sheets <- c('-0.04', '-0.09', '-0.14', '-0.19', '-0.24')
# ---------------------------------------------------------------------------

filename <- file.path(dir.input, "MERFISH_Animal1_cnts.xlsx")
infoname <- file.path(dir.input, "MERFISH_Animal1_info.xlsx")

# MERFISH settings, matching DOST_mHypothalamus.R
neighborhood_threshold <- 3
embedding_dim <- 10

# Number of annotated cell classes
n_cell_classes <- 15

dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

for (sheet in sheets) {

  cnts <- as.data.frame(readxl::read_excel(filename, sheet = sheet))
  row.names(cnts) <- cnts[, "...1"]
  cnts <- cnts[-c(1)]

  xys <- as.data.frame(readxl::read_excel(infoname, sheet = sheet))
  row.names(xys) <- xys[, "...1"]
  gtlabels <- xys$z
  celllabels <- xys$Cell_class
  xys <- xys[-c(1)]
  xys <- xys[-c(-2:-1)]
  xys <- xys[, c(2, 1)]

  sweeps <- list(
    cell_type      = list(R = n_cell_classes, gt = celllabels),
    spatial_domain = list(R = length(unique(gtlabels)), gt = gtlabels)
  )

  for (sweep_name in names(sweeps)) {
    sweep <- sweeps[[sweep_name]]

    results <- list()
    aris <- c()

    for (lambda in lambdas) {
      set.seed(1999)
      res <- DOST(cnts, xys, sweep$R,
                  neighborhood_threshold = neighborhood_threshold,
                  embedding_dim = embedding_dim,
                  lambda = lambda,
                  refinement = FALSE)
      ari <- mclust::adjustedRandIndex(res$labels, sweep$gt)

      results[[as.character(lambda)]] <- res
      aris[as.character(lambda)] <- ari

      cat("section:", sheet, " ", sweep_name, " lambda:", lambda,
          " ARI:", ari, "\n")
    }

    save(results, aris,
         file = file.path(dir.output,
                          paste0("mHypothalamus_DOST_ablation_", sweep_name,
                                 "_sheet_", sheet, ".RData")))
  }
}
