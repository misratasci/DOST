"""
Runtime measurement for STAGATE on DLPFC12.

Counterpart to ADEPT_runtime.py, GraphST_runtime.py, stCluster_runtime.py and
BASS_runtime.R. Times the full pipeline per slice: HVG selection, normalization,
spatial-graph construction, the 1000 training epochs, the neighbors/UMAP step and
mclust. Data loading sits outside the timer.

Adapted from https://github.com/QIFEIDKN/STAGATE_pyG (run_STAGATE.ipynb, runtime cell)
and https://benchmarkst-reproducibility.readthedocs.io/

This is the notebook's runtime cell almost verbatim: it already reloaded the slice
before every repeat and started the timer after the load, so no correction was needed.

Outputs written to dir_output:
  stagate_runtime_<slice>.txt   one wall-clock time in seconds per line, appended
"""

import os
import time
import warnings

import scanpy as sc
import torch

import STAGATE_pyG as STAGATE
from st_loading_utils import load_DLPFC

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Slice indices to time (1 to 12) and how many repeats per slice
slice_indices = range(1, 13)
n_repeats = 5
# ---------------------------------------------------------------------------

#We measure runtime on CPU for fair comparison with other methods
device_name = "cpu"

# Slice IDs
slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

rad_cutoff = 150

device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

for slice_index in slice_indices:
    section_id = str(slices[slice_index - 1])
    n_clusters = 5 if 5 <= slice_index <= 8 else 7

    for repeat in range(n_repeats):
        # Reloaded every repeat, as in the notebook's runtime cell
        ad = load_DLPFC(root_dir=dir_input, section_id=section_id)

        start = time.perf_counter()

        sc.pp.highly_variable_genes(ad, flavor="seurat_v3", n_top_genes=3000)
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        STAGATE.Cal_Spatial_Net(ad, rad_cutoff=rad_cutoff)
        ad = STAGATE.train_STAGATE(ad, device=device)
        sc.pp.neighbors(ad, use_rep="STAGATE")
        sc.tl.umap(ad)
        ad = STAGATE.mclust_R(ad, used_obsm="STAGATE", num_cluster=n_clusters)

        elapsed = time.perf_counter() - start

        print("slice {} repeat {}: {:.2f} s".format(section_id, repeat + 1, elapsed))
        with open(os.path.join(dir_output,
                               "stagate_runtime_{}.txt".format(section_id)),
                  "a+") as fp:
            fp.write(str(elapsed) + "\n")
