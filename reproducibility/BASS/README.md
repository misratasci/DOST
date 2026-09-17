------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# BASS reproducibility

BASS (Li & Zhou, Genome Biology 2022) is a Bayesian hierarchical model that infers cell types and spatial domains jointly: a Potts prior couples neighbouring spots so that domain labels stay spatially smooth, while a mixture model over the expression profiles assigns cell types underneath. It is the only method in this set that returns both, so every script saves `zlabels` (domains) and `clabels` (cell types).

- Original code: <https://github.com/zhengli09/BASS>
- Benchmark protocol: <https://benchmarkst-reproducibility.readthedocs.io/en/latest/BASS_clustering.html>

## Contents

```         
BASS_DLPFC.R            run BASS on DLPFC slices
BASS_BC.R               run BASS on the BC section
BASS_mMAMP.R            run BASS on the mMAMP (MA) section
BASS_mHypothalamus.R    run BASS on the five mHypothalamus sections
BASS_runtime.R          time the full DLPFC pipeline, 5 repeats per slice
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
|------------------------|------------------------|------------------------|
| `BASS_DLPFC.R` | `DLPFC_<slice>_BASS_labels.RData` | `zlabels`, `clabels` |
|  | `bass_aris.RData` | `bass_aris`, one ARI per slice in the loop |
| `BASS_BC.R` | `BC_BASS_labels.RData` | `zlabels`, `clabels`, `ari` |
| `BASS_mMAMP.R` | `mMAMP_BASS_labels.RData` | `zlabels`, `clabels`, `ari` |
| `BASS_mHypothalamus.R` | `BASS_mHypothalamus.RData` | `results` and `aris`, keyed by section |
| `BASS_runtime.R` | `BASS_runtime<index>.txt` | one elapsed time per line, appended |
