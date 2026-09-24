# DOST reproducibility

## Contents

```
DOST_DLPFC.R                        run DOST on DLPFC slices
DOST_BC.R                           run DOST on the BC section
DOST_mMAMP.R                        run DOST on the mMAMP (MA) section
DOST_mHypothalamus.R                run DOST on the five mHypothalamus sections
DOST_runtime.R                      time the DOST call, 5 repeats per slice
DOST_DLPFC_ablation_nGenes.R        sweep the number of HVGs and SVGs
DOST_DLPFC_ablation_embdim.R        sweep the embedding dimension
DOST_DLPFC_ablation_init.R          compare random, PCA and MDS initialization
DOST_DLPFC_ablation_max_iter.R      sweep the number of maximum iterations
DOST_DLPFC_ablation_lambda.R        sweep lambda near the default
DOST_DLPFC_ablation_lambda_sensitivity.R  sweep lambda over 0-1
DOST_mHypothalamus_ablation_lambda.R  sweep lambda, for cell types and domains
load_data.R                         the three Visium dataset loaders
install.R                           installs the R packages the scripts need
README.md
```

## Data

We used the data links provided by the Benchmark ST study: <https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html>

| Dataset | Zenodo | Sections | Domains (R) |
|------------------|------------------|------------------|------------------|
| DLPFC12 | [10698880](https://zenodo.org/records/10698880) | 12 | 7 / 5 / 7 |
| BC | [10698903](https://zenodo.org/records/10698903) | `section1` | 20 |
| mMAMP | [10698931](https://zenodo.org/records/10698931) | `MA` | 52 |
| mHypothalamus | [10698909](https://zenodo.org/records/10698909) | 5 | 8 |

Layout the loaders expect after unzipping:

```
DLPFC12/151673/{spatial/, 151673_filtered_feature_bc_matrix.h5,
                gt/tissue_positions_list_GTs.txt, gt/layered/151673_L1_barcodes.txt, ...}
BC/section1/{spatial/, section1_filtered_feature_bc_matrix.h5,
             gt/gold_metadata.tsv, gt/tissue_positions_list_GTs.txt}
mMAMP/MA/{spatial/, MA_filtered_feature_bc_matrix.h5, metadata.tsv,
          gt/tissue_positions_list_GTs.txt}
mHypothalamus/{MERFISH_Animal1_cnts.xlsx, MERFISH_Animal1_info.xlsx}
```

## Environment

Run `install.R` to install the packages the scripts need.

## Settings

The defaults are tuned for 10x Visium. The non-Visium scripts override them:

| Dataset | `embedding_dim` | `neighborhood_threshold` | `lambda` | `refinement` |
| --- | --- | --- | --- | --- |
| DLPFC12 | 20 | 1 | 0.03 | `TRUE`|
| BC | 20 | 1 | 0.03 | `FALSE` |
| mMAMP | 20 | 1 | 0.03 | `FALSE` |
| mHypothalamus | 10 | 3 | 0.2 | `FALSE` |

`DOST_DLPFC_ablation_init.R` is the one script that does not call `DOST()`:
the exported function always initializes from classical MDS, so the script runs
the same pipeline through the package internals (`preprocess`, `optimize`,
`cluster_embedding`, `refine_labels`) and swaps only the starting point.

### Outputs

| Script | File | Contents |
| --- | --- | --- |
| `DOST_DLPFC.R` | `DOST_DLPFC.RData` | `results` (one entry per slice, each with `Z`, `losses`, `labels`), `aris` |
| | `DOST_DLPFC_noref.RData` | the same, with `refinement = FALSE` |
| `DOST_BC.R` | `BC_DOST_labels.RData` | `labels`, `ari` |
| | `BC_DOST_embeddings.RData` | `emb`, the DOST embedding |
| `DOST_mMAMP.R` | `mMAMP_DOST_labels.RData` | `labels`, `ari` |
| | `mMAMP_DOST_embeddings.RData` | `emb`, the DOST embedding |
| `DOST_mHypothalamus.R` | `mHypothalamus_DOST_sheet_<sheet>_labels.RData` | `labels` |
| | `mHypothalamus_DOST_sheet_<sheet>_embeddings.RData` | `emb` |
| | `mHypothalamus_DOST_aris.RData` | `dost_aris`, one ARI per sheet |
| | `mHypothalamus_DOST_celltype_sheet_<sheet>_labels.RData` | `labels`, the `lambda = 0` run |
| | `mHypothalamus_DOST_celltype_sheet_<sheet>_embeddings.RData` | `emb` |
| | `mHypothalamus_DOST_celltype_aris.RData` | `dost_celltype_aris` |
| `DOST_runtime.R` | `DOST_runtime_slice<index>.txt` | one elapsed time per line, appended |
| `DOST_DLPFC_ablation_nGenes.R` | `DOST_DLPFC_ablation_nGenes_HVG.RData` | `aris` (slices x gene counts), `ari_means` |
| | `DOST_DLPFC_ablation_nGenes_SVG.RData` | the same, with SVG selection |
| `DOST_DLPFC_ablation_embdim.R` | `DOST_DLPFC<index>_embdim.RData` | `results` keyed by embedding dim, `aris` |
| `DOST_DLPFC_ablation_init.R` | `DOST_DLPFC<index>_init_<init>.RData` | `results`, `ari` |
| `DOST_DLPFC_ablation_max_iter.R` | `DOST_DLPFC_ablation_max_iter_<index>.RData` | `results` keyed by iteration budget |
| `DOST_DLPFC_ablation_lambda.R` | `DOST_DLPFC_ablation_lambda_<index>.RData` | `results` keyed by lambda |
| `DOST_DLPFC_ablation_lambda_sensitivity.R` | `DOST_DLPFC<index>_lambda_noref.RData` | `results` keyed by lambda, `aris` |
| | `DOST_DLPFC<index>_lambda_ref.RData` | the same, labels refined |
| `DOST_mHypothalamus_ablation_lambda.R` | `mHypothalamus_DOST_ablation_cell_type_sheet_<sheet>.RData` | `results` keyed by lambda, `aris` |
| | `mHypothalamus_DOST_ablation_spatial_domain_sheet_<sheet>.RData` | the same, for domains |
