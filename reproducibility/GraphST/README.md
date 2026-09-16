# GraphST reproducibility

GraphST (Long *et al.*, Nature Communications 2023) is a graph self-supervised
contrastive learning model: it learns spot embeddings from expression and spatial
position, clusters them with `mclust`, then refines each label by majority vote over
the spot's nearest neighbours.

- Original code: <https://github.com/JinmiaoChenLab/GraphST>
- Benchmark protocol: <https://benchmarkst-reproducibility.readthedocs.io/>

## Contents

```
GraphST_DLPFC.py            run GraphST on one DLPFC slice
GraphST_BC.py               run GraphST on the BC section
GraphST_mMAMP.py            run GraphST on the mMAMP (MA) section
GraphST_mHypothalamus.py    run GraphST on one mHypothalamus section
GraphST_runtime.py          time the full DLPFC pipeline, 5 repeats per slice
st_loading_utils.py         dataset loaders of BenchmarkST, unmodified
GraphST/                    the GraphST package, unmodified from the original repository
environment.yml             conda environment
```

## Data

We used the data links provided by the BenchmarkST study:
<https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html>

| Dataset | Zenodo | Sections | Domains |
| --- | --- | --- | --- |
| DLPFC12 | [10698880](https://zenodo.org/records/10698880) | 12 | 7 / 5 / 7 |
| BC | [10698903](https://zenodo.org/records/10698903) | `section1` | 20 |
| mMAMP | [10698931](https://zenodo.org/records/10698931) | `MA` | 52 |
| mHypothalamus | [10698909](https://zenodo.org/records/10698909) | 5 | 8 |

Layout the loaders expect after unzipping:

```
DLPFC12/151673/{spatial/, 151673_filtered_feature_bc_matrix.h5, gt/tissue_positions_list_GTs.txt}
BC/section1/{spatial/, section1_filtered_feature_bc_matrix.h5, gt/tissue_positions_list_GTs.txt}
mMAMP/MA/{spatial/, MA_filtered_feature_bc_matrix.h5, gt/tissue_positions_list_GTs.txt}
mHypothalamus/{MERFISH_Animal1_cnts.xlsx, MERFISH_Animal1_info.xlsx}
```

## Environment

```
conda env create -f environment.yml
conda activate graphst
```

### Outputs

| Script | Labels | Embedding | PCA embedding | ARI |
| --- | --- | --- | --- | --- |
| `GraphST_DLPFC.py` | `graphst_slice_<slice>.csv` | `graphst_embeddings_slice_<slice>.csv` | `graphst_embeddings_pca_slice_<slice>.csv` | `graphst_aris.txt` |
| `GraphST_BC.py` | `graphst_bc.csv` | `graphst_embeddings_bc.csv` | `graphst_embeddings_pca_bc.csv` | `graphst_bc_aris.txt` |
| `GraphST_mMAMP.py` | `graphst_mMAMP.csv` | `graphst_embeddings_mMAMP.csv` | `graphst_embeddings_pca_mMAMP.csv` | `graphst_mMAMP_aris.txt` |
| `GraphST_mHypothalamus.py` | `graphst_mHypothalamus_<section>.csv` | `graphst_embeddings_mHypothalamus_<section>.csv` | `graphst_embeddings_pca_mHypothalamus_<section>.csv` | `graphst_aris_mHypothalamus.txt` |

`GraphST_runtime.py` writes only `graphst_runtime_<slice>.txt`, one wall-clock time in
seconds per line, appended across repeats.