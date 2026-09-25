# Runtime figures for DLPFC12.
#
# Reads the files the scripts in reproducibility/<method>/ write, so run those
# first with their dir.output all pointing at the same folder, then set
# dir.output below to that folder.
#
# Inputs (all from dir.output). The R methods index by slice position 1-12,
# the Python methods by slice ID:
#   BASS       BASS_runtime<index>.txt
#   ADEPT      adept_runtime_<slice>.txt
#   GraphST    graphst_runtime_<slice>.txt
#   stCluster  stcluster_runtime_<slice>.txt
#   STAGATE    stagate_runtime_<slice>.txt
#   SpaGCN     spagcn_runtime_<slice>.txt
#   BANKSY     BANKSY_runtime<index>.txt
#   DR.SC      DRSC_runtime<index>.txt
#   DOST       DOST_runtime_slice<index>.txt
#
# Outputs written to dir.figures:
#   DLPFC_runtime_violinplot.pdf         mean runtime per slice, one violin per method
#   supplementary_all_slices_runtime.pdf mean +/- sd per method, faceted by slice
#   runtime_detailed_table.csv           the same numbers as "mean +/- sd"

source("reproducibility/figures/figure_utils.R")

library(dplyr)

# ---------------------------------------------------------------------------
# Change this path to the folder the runtime scripts wrote their results to
dir.output <- "path/to/output"

# Change this path to where you want the figures saved
dir.figures <- "path/to/figures"
# ---------------------------------------------------------------------------

# Slice IDs
slices <- c(151507:151510, 151669:151676)

# Order the methods appear in, in every runtime panel
method_order <- c("BASS", "ADEPT", "GraphST", "stCluster",
                  "STAGATE", "SpaGCN", "BANKSY", "DR.SC", "DOST")

dir.create(dir.figures, showWarnings = FALSE, recursive = TRUE)

runtime_paths <- function(slice_index) {
  slice_id <- slices[slice_index]
  list(
    BASS      = file.path(dir.output, paste0("BASS_runtime", slice_index, ".txt")),
    ADEPT     = file.path(dir.output, paste0("adept_runtime_", slice_id, ".txt")),
    GraphST   = file.path(dir.output, paste0("graphst_runtime_", slice_id, ".txt")),
    stCluster = file.path(dir.output, paste0("stcluster_runtime_", slice_id, ".txt")),
    STAGATE   = file.path(dir.output, paste0("stagate_runtime_", slice_id, ".txt")),
    SpaGCN    = file.path(dir.output, paste0("spagcn_runtime_", slice_id, ".txt")),
    BANKSY    = file.path(dir.output, paste0("BANKSY_runtime", slice_index, ".txt")),
    DR.SC     = file.path(dir.output, paste0("DRSC_runtime", slice_index, ".txt")),
    DOST      = file.path(dir.output, paste0("DOST_runtime_slice", slice_index, ".txt"))
  )
}

load_all_runs <- function(slice_index) {
  paths <- runtime_paths(slice_index)
  rows <- lapply(names(paths), function(m) {
    if (!file.exists(paths[[m]])) return(NULL)
    vals <- read.csv(paths[[m]], header = FALSE)$V1
    data.frame(method = m, runtime = vals, slice_id = slices[slice_index])
  })
  bind_rows(rows)
}

all_data <- bind_rows(lapply(seq_along(slices), load_all_runs))
if (nrow(all_data) == 0) {
  stop("no runtime files found in ", dir.output)
}
all_data$method <- factor(all_data$method, levels = method_order)

missing <- setdiff(method_order, unique(as.character(all_data$method)))
if (length(missing) > 0) {
  message("no runtime files for: ", paste(missing, collapse = ", "))
}

summary_per_slice <- all_data %>%
  group_by(slice_id, method) %>%
  summarise(mean_runtime = mean(runtime),
            sd_runtime = sd(runtime),
            .groups = "drop")

# ---------------------------------------------------------------------------
# Main figure: distribution of the per-slice means, log scale
# ---------------------------------------------------------------------------

violin_plot <- ggplot(summary_per_slice,
                      aes(x = method, y = mean_runtime, fill = method)) +
  geom_violin(trim = FALSE, color = "black", width = 0.8, alpha = 1) +
  ggbeeswarm::geom_beeswarm(cex = 1.0, size = 0.5, priority = "density",
                            alpha = 1, color = "black") +
  scale_y_log10() +
  guides(fill = "none", color = "none") +
  labs(y = "Mean Runtime (s)", x = NULL) +
  theme_classic()

save_cropped(file.path(dir.figures, "DLPFC_runtime_violinplot.pdf"),
             violin_plot, width = 6, height = 4)

# ---------------------------------------------------------------------------
# Supplementary figure: per-slice bars with the spread over the repeats
# ---------------------------------------------------------------------------

supp_plot <- ggplot(summary_per_slice,
                    aes(x = method, y = mean_runtime, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_runtime - sd_runtime,
                    ymax = mean_runtime + sd_runtime),
                width = 0.3, linewidth = 0.6) +
  # An explicit line at zero so every facet has an x-axis, not just the bottom row
  geom_hline(yintercept = 0, linewidth = 0.5) +
  facet_wrap(~slice_id, ncol = 4) +
  labs(x = NULL, y = "Mean Runtime (s)") +
  guides(fill = "none") +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_blank())

save_cropped(file.path(dir.figures, "supplementary_all_slices_runtime.pdf"),
             supp_plot, width = 12, height = 9)

# ---------------------------------------------------------------------------
# Supplementary table
# ---------------------------------------------------------------------------

summary_stats <- summary_per_slice %>%
  group_by(method) %>%
  summarise(avg_runtime_sec = mean(mean_runtime),
            sd_runtime_sec = sd(mean_runtime),
            min_runtime = min(mean_runtime),
            max_runtime = max(mean_runtime)) %>%
  arrange(factor(method, levels = method_order))
print(summary_stats)

# Slices as rows, methods as columns, "mean +/- sd" over the repeats
runtime_table <- summary_per_slice %>%
  mutate(val = paste0(round(mean_runtime, 2), " ± ", round(sd_runtime, 2))) %>%
  select(slice_id, method, val) %>%
  tidyr::pivot_wider(names_from = method, values_from = val)

write.csv(runtime_table, file.path(dir.figures, "runtime_detailed_table.csv"),
          row.names = FALSE)
print(runtime_table)

# LaTeX for the supplementary table
print(knitr::kable(runtime_table, format = "latex", booktabs = TRUE,
                   caption = "Runtime Details"))
