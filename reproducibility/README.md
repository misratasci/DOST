# Reproducibility

Scripts that reproduce the results and figures of the paper.

```
DOST/        DOST, its ablation studies, and the runtime benchmark (R)
BANKSY/      baseline, R
BASS/        baseline, R
DR.SC/       baseline, R
GraphST/     baseline, Python (clone of the original code needed, see its README)
STAGATE/     baseline, Python
SpaGCN/      baseline, Python
stCluster/   baseline, Python
ADEPT/       baseline, Python (clone of the original code plus a patch needed, see its README)
figures/     draws the figures from the files the scripts above write
licenses/    license texts of third-party code
```

## Order

1. Download the datasets, as described in `DOST/README.md`.
2. Run the method scripts. Each folder's README lists its scripts, environment and output
   files, and the scripts write to the path you set at the top of each file.
3. Run the scripts in `figures/`, which read those outputs.

## Third-party code

Parts of this directory were written by others and keep their own licenses. See
`../THIRD_PARTY_NOTICES.md`.
