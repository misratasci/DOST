# SpaGCN reproducibility

SpaGCN (Hu *et al.*, Nature Methods 2021) builds a weighted graph over spots from their
spatial coordinates — optionally weighted by local histology colour — runs a graph
convolutional autoencoder with a deep embedded clustering head initialised by louvain,
and then refines each label by majority vote over its nearest neighbours.

- Original code: <https://github.com/jianhuupenn/SpaGCN>
- Benchmark protocol: <https://benchmarkst-reproducibility.readthedocs.io/>

## Contents

```
SpaGCN_DLPFC.py            run SpaGCN on one DLPFC slice
SpaGCN_BC.py               run SpaGCN on the BC section
SpaGCN_mMAMP.py            run SpaGCN on the mMAMP (MA) section
SpaGCN_mHypothalamus.py    run SpaGCN on one mHypothalamus section
SpaGCN_runtime.py          time the full DLPFC pipeline, 5 repeats per slice
st_loading_utils.py        dataset loaders from BenchmarkST, unmodified
SpaGCN/                    the SpaGCN package, as used for these runs
environment.yml            conda environment
```

`SpaGCN/` carries a **deliberate change** to `SpaGCN.py`, 
described under Notes below: `SpaGCN.predict` optionally takes and
returns an AnnData so the latent embedding can be saved.

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
conda activate spagcn
```

### Outputs

| Script | Labels | Embedding | ARI |
| --- | --- | --- | --- |
| `SpaGCN_DLPFC.py` | `spagcn_slice_<slice>.csv` | `spagcn_embeddings_slice_<slice>.csv` | `spagcn_aris.txt` |
| `SpaGCN_BC.py` | `spagcn_bc.csv` | `spagcn_embeddings_bc.csv` | `spagcn_aris_bc.txt` |
| `SpaGCN_mMAMP.py` | `spagcn_mMAMP.csv` | `spagcn_embeddings_mMAMP.csv` | `spagcn_aris_mMAMP.txt` |
| `SpaGCN_mHypothalamus.py` | `spagcn_mHypothalamus_<section>.csv` | `spagcn_embeddings_mHypothalamus_<section>.csv` | `spagcn_aris_mHypothalamus.txt` |

`SpaGCN_runtime.py` writes only `spagcn_runtime_<slice>.txt`, one wall-clock time in
seconds per line, appended across repeats.

## Notes

**`SpaGCN.predict` was modified to expose the embedding.** Upstream it is
`predict(self)`, returning only the cluster assignments and the per-class probabilities
and throwing away `z`, the latent matrix the model computes — which is why stock SpaGCN
gives you labels but no embedding. Here it is `predict(self, adata=None)`: called with an
AnnData it stores `z` as `adata.obsm['X_spagcn']` and returns `(y_pred, prob, adata)`;
called with no argument it behaves exactly as before and returns `(y_pred, prob)`. The
argument has to stay optional, because `util.search_res` calls `clf.predict()` with no
arguments while searching for the resolution, so both call sites must keep working.
