# STAGATE reproducibility

STAGATE (Dong & Zhang, Nature Communications 2022) is a graph attention autoencoder:
it builds a spatial neighbour graph from a distance cut-off, learns a low-dimensional
embedding by reconstructing the expression matrix through attention-weighted
neighbourhoods, and clusters that embedding with `mclust`.

- Original code: <https://github.com/QIFEIDKN/STAGATE_pyG>
- Benchmark protocol: <https://benchmarkst-reproducibility.readthedocs.io/>

## Contents

```
STAGATE_DLPFC.py            run STAGATE on one DLPFC slice
STAGATE_BC.py               run STAGATE on the BC section
STAGATE_mMAMP.py            run STAGATE on the mMAMP (MA) section
STAGATE_mHypothalamus.py    run STAGATE on one mHypothalamus section
STAGATE_runtime.py          time the full DLPFC pipeline, 5 repeats per slice
st_loading_utils.py         the benchmark's dataset loaders, unmodified
STAGATE_pyG/                the STAGATE package, cloned from original repository, unmodified
environment.yml             conda environment
```

## Data

We used the data links provided by the Benchmark ST study:
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
conda activate stagate
```

### Outputs

| Script | Labels | Embedding | ARI |
| --- | --- | --- | --- |
| `STAGATE_DLPFC.py` | `stagate_slice_<slice>.csv` | `stagate_embeddings_slice_<slice>.csv` | `stagate_aris.txt` |
| `STAGATE_BC.py` | `stagate_bc.csv` | `stagate_embeddings_bc.csv` | `stagate_aris_bc.txt` |
| `STAGATE_mMAMP.py` | `stagate_mMAMP.csv` | `stagate_embeddings_mMAMP.csv` | `stagate_aris_mMAMP.txt` |
| `STAGATE_mHypothalamus.py` | `stagate_mHypothalamus_<section>.csv` | `stagate_embeddings_mHypothalamus_<section>.csv` | `stagate_aris_mHypothalamus.txt` |

`STAGATE_runtime.py` writes only `stagate_runtime_<slice>.txt`, one wall-clock time in
seconds per line, appended across repeats.
