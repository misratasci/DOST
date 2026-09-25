# BANKSY reproducibility

BANKSY (Singhal *et al.*, Nature Genetics 2024) augments each cell's expression vector
with two neighbourhood summaries — the local mean and an azimuthal Gabor filter — then
runs PCA and leiden clustering on the augmented matrix. A single parameter, `lambda`,
sets how much weight the neighbourhood terms carry, which is what lets the same method
serve both cell typing (low lambda) and domain segmentation (higher lambda).

- Original code: <https://prabhakarlab.github.io/Banksy/>
- Preprocessing follows <https://prabhakarlab.github.io/Banksy/articles/multi-sample.html>

## Contents

```
BANKSY_DLPFC.R            run BANKSY on DLPFC slices
BANKSY_BC.R               run BANKSY on the BC section
BANKSY_mMAMP.R            run BANKSY on the mMAMP (MA) section
BANKSY_mHypothalamus.R    run BANKSY on the five mHypothalamus sections
BANKSY_runtime.R          time the BANKSY pipeline, 5 repeats per slice
load_data.R               the three dataset loaders
install.R                 installs the R packages the scripts need
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

Run install.R to install the packages the scripts need.

### Outputs

| Script | File | Contents |
| --- | --- | --- |
| `BANKSY_DLPFC.R` | `DLPFC_BANKSY_slice_<slice>_labels.RData` | `labels` |
| | `DLPFC_BANKSY_slice_<slice>_embeddings.RData` | `emb`, the BANKSY embedding |
| | `DLPFC_BANKSY_aris.RData` | `banksy_aris`, one ARI per slice in the loop |
| `BANKSY_BC.R` | `BC_BANKSY_labels.RData` | `labels` |
| | `BC_BANKSY_embeddings.RData` | `emb`, the the BANKSY embedding |
| `BANKSY_mMAMP.R` | `mMAMP_BANKSY_labels.RData` | `labels` |
| | `mMAMP_BANKSY_embeddings.RData` | `emb`, the the BANKSY embedding |
| `BANKSY_mHypothalamus.R` | `mHypothalamus_BANKSY_sheet_<sheet>_labels.RData` | `labels` |
| | `mHypothalamus_BANKSY_sheet_<sheet>_embeddings.RData` | `emb`, the BANKSY embedding |
| | `mHypothalamus_BANKSY_aris.RData` | `banksy_aris`, one ARI per sheet in the loop |
| `BANKSY_runtime.R` | `BANKSY_runtime<index>.txt` | one elapsed time per line, appended |

