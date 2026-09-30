# BC (human breast cancer, Visium section1) lambda sensitivity figure for DOST.
#
# Reads the file DOST_BC_ablation_lambda_sensitivity.R writes in
# reproducibility/DOST/, so run that first, then set dir.output below to the
# folder it wrote to.
#
# Inputs (all from dir.output):
#   DOST_BC_ablation_lambda_sensitivity.R
#                                 DOST_BC_lambda_sensitivity.RData     results, aris
#
# Outputs written to dir.figures:
#   BC_lambda_sensitivity.pdf      ARI against lambda
#   BC_lambda_panel.pdf            clustering and UMAP per lambda

source("reproducibility/figures/load_data.R")
source("reproducibility/figures/figure_utils.R")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir.input <- "path/to/data/BC"

# Change this path to the folder the ablation script wrote its results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <- "path/to/figures"
# ---------------------------------------------------------------------------

set.seed(1999)

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

results_file <- file.path(dir.output, "DOST_BC_lambda_sensitivity.RData")
results <- load_object(results_file, "results")
aris <- load_object(results_file, "aris")

# ---------------------------------------------------------------------------
# Lambda sensitivity
# ---------------------------------------------------------------------------

lambda_results <- data.frame(Lambda = as.numeric(names(aris)),
                             ARI = as.numeric(aris))

lambda_plot <- ggplot(lambda_results, aes(x = Lambda, y = ARI)) +
  geom_line(linewidth = 0.8, color = "#E41A1C") +
  geom_point(size = 2, color = "#E41A1C") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_x_continuous(breaks = lambda_results$Lambda) +
  labs(x = "Lambda (Spatial Regularization)",
       y = "Adjusted Rand Index (ARI)") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_blank())

save_cropped(file.path(dir.figures, "BC_lambda_sensitivity.pdf"),
             lambda_plot, width = 5, height = 3.5)
