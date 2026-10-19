# Third-party notices

The DOST R package (`R/`, `src/`, `man/`, `vignettes/`) and the scripts
written for this project are released under GPL-3 (see `LICENSE.md`).

The `reproducibility/` directory also contains code written by others,
used to run the baseline methods we compare against. That code keeps its
original license and copyright. Each component is listed below.

## Included under the MIT License

The full license text, with the original copyright line, is in the
folder named below.

| Component | Location | Source | Copyright | Changes |
|----|----|----|----|----|
| SpaGCN | `reproducibility/SpaGCN/SpaGCN/` | <https://github.com/jianhuupenn/SpaGCN> | \(c\) 2020 JianHu | `SpaGCN.predict` extended to return the embedding, see `reproducibility/SpaGCN/README.md` |
| STAGATE (PyG) | `reproducibility/STAGATE/STAGATE_pyG/` | <https://github.com/QIFEIDKN/STAGATE_pyG> (now <https://github.com/RucDongLab/STAGATE_pyG>) | \(c\) 2021 Kangning Dong | None |
| stCluster | `reproducibility/stCluster/stCluster/` | <https://github.com/hannshu/stCluster> | \(c\) 2023 han shu | See `reproducibility/stCluster/README.md` |
| BenchmarkST dataset loaders (`st_loading_utils.py`) | `reproducibility/{GraphST,STAGATE,SpaGCN,stCluster}/st_loading_utils.py` | <https://github.com/maiziezhoulab/BenchmarkST> | \(c\) 2021 MaizieZhouLab | Adapted in `stCluster`, unmodified elsewhere |

The BenchmarkST license text is in
`reproducibility/licenses/BenchmarkST_LICENSE.txt`.

## Not redistributed

The following baseline methods are not part of this archive. Obtain them
from their original repositories; the reproducibility READMEs describe
how they are used.

| Component | Source                                      | License             |
|-----------|---------------------------------------------|---------------------|
| GraphST   | <https://github.com/JinmiaoChenLab/GraphST> | AGPL-3.0            |
| ADEPT     | <https://github.com/maiziezhoulab/ADEPT>    | No license declared |

`reproducibility/ADEPT/ADEPT_compat.patch` is our own small patch to
ADEPT, to be applied to a clone of the original repository.

## R baselines

BANKSY, BASS and DR.SC are installed from CRAN, Bioconductor or GitHub
by the `install.R` scripts and are not redistributed here. The same
applies to every package listed under `Imports` and `Suggests` in
`DESCRIPTION`.

## Data

No datasets are included. The datasets used in the paper are obtained
from the sources listed in `reproducibility/DOST/README.md`, under their
own terms.
