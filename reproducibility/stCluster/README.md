# stCluster reproducibility

stCluster (Shu *et al.*, 2024) is a graph contrastive learning model: it builds a
spatial graph, prunes edges between spots that a louvain pre-clustering puts in
different domains, and trains a graph attention autoencoder under a combined
contrastive, reconstruction, adjacency and DEC-prediction loss. The embedding is then
clustered with `mclust`.

- Original code: <https://github.com/hannshu/stCluster>

## Contents

```
stCluster_DLPFC.py          run stCluster on one DLPFC slice
stCluster_BC.py             run stCluster on the BC section
stCluster_mMAMP.py          run stCluster on the mMAMP (MA) section
stCluster_mHypothalamus.py  run stCluster on one mHypothalamus section
stCluster_runtime.py        time the full DLPFC pipeline, 5 repeats per slice
st_loading_utils.py         the four dataset loaders, adapted from BenchmarkST (Hu et al.)
stCluster/                  the stCluster package
environment.yml             conda environment
```

`run.py`'s `mclust_R` in `stCluster/` was rewritten locally to
work around an rpy2 argument-wrapping error, and the original implementation is still
in the file as unreachable code after the `return`.

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
conda activate stcluster
```

### Outputs

| Script | Labels | Embedding | ARI |
| --- | --- | --- | --- |
| `stCluster_DLPFC.py` | `stcluster_slice_<slice>.csv` | `stcluster_embeddings_slice_<slice>.csv` | `stcluster_aris.txt` |
| `stCluster_BC.py` | `stcluster_bc_section1.csv` | `stcluster_embeddings_bc_section1.csv` | `stcluster_aris_bc.txt` |
| `stCluster_mMAMP.py` | `stcluster_mmamp_MA.csv` | `stcluster_embeddings_mmamp_MA.csv` | `stcluster_aris_mmamp.txt` |
| `stCluster_mHypothalamus.py` | `stcluster_mhypo_<section>.csv` | `stcluster_embeddings_mhypo_<section>.csv` | `stcluster_aris_mhypo.txt` |
