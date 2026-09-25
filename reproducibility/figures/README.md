# Article figures

These scripts draw the figures in the article. They do not run any method: they
read the files the scripts in `reproducibility/<method>/` write and turn them
into panels.

## Contents

```
DLPFC_figure.R                  slice grid and ARI violin, DLPFC12
DLPFC_ablation_figure.R         gene count, embedding dim, init, max iter, lambda
BC_figure.R                     spatial and UMAP rows, BC section1
mMAMP_figure.R                  spatial and UMAP rows, mMAMP MA
mHypothalamus_figure.R          section grid and ARI plots, mHypothalamus
mHypothalamus_ablation_figure.R lambda sweep, cell types against domains
runtime_figure.R                runtime comparison and supplementary table
figure_utils.R                  readers for the method outputs, panel styles
load_data.R                     dataset loaders
install.R                       installs the R packages the figure scripts need
README.md
```

## Environment

Run `install.R`. `pdfcrop` comes from TeX Live or MacTeX and is optional.

## Figures

| Script | File | Contents |
| --- | --- | --- |
| `DLPFC_figure.R` | `DLPFC_slices<selected>.pdf` | ground truth and all ten columns, selected slices |
| | `DLPFC_violinplot.pdf` | ARI over the 12 slices, one violin per method |
| `DLPFC_ablation_figure.R` | `nGenes.pdf` | mean ARI against HVG and SVG count |
| | `DLPFC_ablation_embedding_dim_<index>.pdf` | clustering and UMAP per embedding dim |
| | `DLPFC_initialization_<index>.pdf` | clustering and UMAP per initialization |
| | `DLPFC_ablation_max_iter_<index>.pdf` | clustering and UMAP per iteration budget |
| | `DLPFC_ablation_lambda_<index>.pdf` | clustering and UMAP per lambda |
| | `DLPFC_lambda_sensitivity.pdf` | ARI against lambda, faceted by slice |
| `BC_figure.R` | `BC.pdf` | spatial row and UMAP row |
| `mMAMP_figure.R` | `mMAMP.pdf` | spatial row and UMAP row |
| `mHypothalamus_figure.R` | `mHypo_slices.pdf` | ground truth and all methods, five sections |
| | `mHypo_violinplot.pdf` | ARI over the five sections |
| `mHypothalamus_ablation_figure.R` | `mHypothalamus_ablation_<sheet>.pdf` | lambda sweep, for one chosen section |
| `runtime_figure.R` | `DLPFC_runtime_violinplot.pdf` | mean runtime per slice, log scale |
| | `supplementary_all_slices_runtime.pdf` | mean +/- sd per method, faceted by slice |
| | `runtime_detailed_table.csv` | the same numbers as "mean +/- sd" |

Every PDF is passed through `pdfcrop` after `ggsave`; without `pdfcrop` on the
path the plots are still written, just uncropped.
