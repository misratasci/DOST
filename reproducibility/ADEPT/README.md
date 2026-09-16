# ADEPT reproducibility

ADEPT (Hu *et al.*, RECOMB-seq 2023) is a graph autoencoder that clusters spatial
transcriptomics data by selecting differentially expressed genes (DEGs) and imputing
the resulting DEG matrices before a final clustering step.

- Original code: <https://github.com/maiziezhoulab/ADEPT>
- Benchmark: <https://benchmarkst-reproducibility.readthedocs.io/en/latest/ADEPT.html/>

## Contents

```
ADEPT_DLPFC.py           run ADEPT on one DLPFC slice
ADEPT_BC.py              run ADEPT on the BC section
ADEPT_mMAMP.py           run ADEPT on the mMAMP (MA) section
ADEPT_mHypothalamus.py   run ADEPT on one mHypothalamus section
ADEPT_runtime.py         time the full DLPFC pipeline, 5 repeats per slice
GAAE/                    the ADEPT package, unmodified from the original repository
environment.yml          conda environment
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

Layout ADEPT's loaders expect after unzipping:

```
DLPFC12/151673/{spatial/, 151673_filtered_feature_bc_matrix.h5, gt/tissue_positions_list_GTs.txt}
BC/section1/{spatial/, section1_filtered_feature_bc_matrix.h5, gt/tissue_positions_list_GTs.txt}
mMAMP/MA/{spatial/, MA_filtered_feature_bc_matrix.h5, gt/tissue_positions_list_GTs.txt}
mHypothalamus/{MERFISH_Animal1_cnts.xlsx, MERFISH_Animal1_info.xlsx}
```

## Environment

```
conda env create -f environment.yml
conda activate adept
```

### Outputs

| Script | Labels | Embedding | ARI |
| --- | --- | --- | --- |
| `ADEPT_DLPFC.py` | `DLPFC_<slice>_adept.txt` | `adept_embeddings_slice_<slice>.csv` | `adept_aris.txt` |
| `ADEPT_BC.py` | `BC_adept.txt` | `adept_embeddings_BC.csv` | `adept_aris_bc.txt` |
| `ADEPT_mMAMP.py` | `mMAMP_adept.txt` | `adept_embeddings_mMAMP.csv` | `adept_aris_mMAMP.txt` |
| `ADEPT_mHypothalamus.py` | `mHypothalamus_<section>_adept.txt` | `adept_embeddings_mHypothalamus_<section>.csv` | `adept_aris_mHypothalamus.txt` |

`ADEPT_runtime.py` writes only `adept_runtime_<slice>.txt`, one wall-clock time in
seconds per line, appended across repeats.