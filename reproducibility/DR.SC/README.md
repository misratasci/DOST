# DR.SC reproducibility

DR.SC (Liu *et al.*, Nucleic Acids Research 2022) fits dimension reduction and spatial
clustering as one model rather than in sequence: a hidden Markov random field over the
spot neighbourhood graph supplies the spatial prior, and the low-dimensional embedding
and the domain labels are estimated jointly by an EM algorithm. Because the two steps
are not separated, the embedding it returns is the one the clustering actually used.

- Original code: <https://feiyoung.github.io/DR.SC/>
- Benchmark protocol: <https://benchmarkst-reproducibility.readthedocs.io/>

## Contents

```
DRSC_DLPFC.R            run DR.SC on DLPFC slices
DRSC_BC.R               run DR.SC on the BC section
DRSC_mMAMP.R            run DR.SC on the mMAMP (MA) section
DRSC_mHypothalamus.R    run DR.SC on the five mHypothalamus sections
DRSC_runtime.R          time the DR.SC pipeline, 5 repeats per slice
load_data.R             the three dataset loaders
install.R               installs the R packages the scripts need
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
| `DRSC_DLPFC.R` | `DLPFC_DRSC_labels_<slice>.RData` | `labels`, the `spatial.drsc.cluster` assignment |
| | `DLPFC_DRSC_embeddings_<slice>.RData` | `emb`, the DR.SC embedding |
| | `DLPFC_DRSC_aris.RData` | `drsc_aris`, one ARI per slice in the loop |
| `DRSC_BC.R` | `BC_DRSC_labels.RData` | `labels` |
| | `BC_DRSC_embeddings.RData` | `emb`, the DR.SC embedding |
| `DRSC_mMAMP.R` | `mMAMP_DRSC_labels.RData` | `labels` |
| | `mMAMP_DRSC_embeddings.RData` | `emb`, the DR.SC embedding |
| `DRSC_mHypothalamus.R` | `mHypothalamus_DRSC_sheet_<sheet>_labels.RData` | `labels` |
| | `mHypothalamus_DRSC_sheet_<sheet>_embeddings.RData` | `emb`, the DR.SC embedding |
| | `mHypothalamus_DRSC_aris.RData` | `drsc_aris`, one ARI per sheet in the loop |
| `DRSC_runtime.R` | `DRSC_runtime<index>.txt` | one elapsed time per line, appended |

