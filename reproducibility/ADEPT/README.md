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
ADEPT_compat.patch       our small compatibility patch for the original code (see Setup)
ADEPT_src/               not included: clone of the original repository (see Setup)
environment.yml          conda environment
```

## Setup

The ADEPT repository has no license, so its code is not redistributed here. Clone the
original repository next to the scripts at the commit we used, then apply our patch:

```
cd reproducibility/ADEPT
git clone https://github.com/maiziezhoulab/ADEPT ADEPT_src
git -C ADEPT_src checkout 149c3e964dd392787c4c9710a7137a1e69d0419f
git -C ADEPT_src apply ../ADEPT_compat.patch
```

The patch changes two things in `GAAE/`: it replaces a `set_axis(..., inplace=False)` call
and a `var.index[:, np.newaxis]` indexing that recent pandas versions reject, and it routes
the `MA` section to the `load_mMAMP` loader. The scripts add `ADEPT_src/` to `sys.path` and
import `GAAE` from there.

## Data

We used the data links provided by the BenchmarkST study:
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

## License

ADEPT has no declared license and is not redistributed here; `ADEPT_compat.patch` is our own change to it. See also `THIRD_PARTY_NOTICES.md` in the repository root.
