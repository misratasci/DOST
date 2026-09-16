"""
Runtime measurement for GraphST on DLPFC12.

Adapted from https://benchmarkst-reproducibility.readthedocs.io/

Times the full GraphST pipeline repeated `n_repeats` times per slice. Nothing is saved except the timings.

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download DLPFC12 data from https://zenodo.org/records/10698880

Outputs written to dir_output:
  graphst_runtime_<slice>.txt   one wall-clock time in seconds per line, appended
"""

import os
import time
import warnings

import torch

from GraphST import GraphST
from GraphST.utils import clustering
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

# We measure runtime with CPU for fair comparison with other CPU-only methods
device_name = "cpu"

slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

radius = 50
tool = "mclust"

device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

for slice_index in slice_indices:
    section_id = str(slices[slice_index - 1])
    n_clusters = 5 if 5 <= slice_index <= 8 else 7

    for repeat in range(n_repeats):
        ad = load_DLPFC(root_dir=dir_input, section_id=section_id)

        start = time.perf_counter()

        model = GraphST.GraphST(ad, device=device)
        ad_out = model.train()
        if tool == "mclust":
            clustering(ad_out, n_clusters, radius=radius, method=tool, refinement=True)
        elif tool in ["leiden", "louvain"]:
            clustering(ad_out, n_clusters, radius=radius, method=tool,
                       start=0.1, end=2.0, increment=0.01, refinement=False)

        elapsed = time.perf_counter() - start

        print("slice {} repeat {}: {:.2f} s".format(section_id, repeat + 1, elapsed))
        with open(os.path.join(dir_output,
                               "graphst_runtime_{}.txt".format(section_id)), "a+") as fp:
            fp.write(str(elapsed) + "\n")
