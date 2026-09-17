"""
Runtime measurement for stCluster on DLPFC12.

Times the full pipeline per slice: preprocessing, spatial-graph construction,
louvain pre-clustering, the training loop and the mclust clustering.
Data loading sits outside the timer.

Adapted from https://github.com/hannshu/stCluster

Outputs written to dir_output:
  stcluster_runtime_<slice>.txt   one wall-clock time in seconds per line, appended
"""

import os
import time
import warnings

import scanpy as sc

from st_loading_utils import load_DLPFC
from stCluster.run import evaluate_embedding
from stCluster.train import train

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

slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

os.makedirs(dir_output, exist_ok=True)

for slice_index in slice_indices:
    section_id = str(slices[slice_index - 1])
    n_cluster = 5 if 5 <= slice_index <= 8 else 7

    for repeat in range(n_repeats):
        ad = load_DLPFC(root_dir=dir_input, section_id=section_id)
        
        start = time.perf_counter()

        sc.pp.highly_variable_genes(ad, flavor="seurat_v3", n_top_genes=3000)
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)

        ad, g = train(ad, radius=150, ae_rate=0.8, adj_rate=0.2, pred_rate=0.3, seed=0)
        ad, score = evaluate_embedding(adata=ad, n_cluster=n_cluster,
                                       cluster_method=["mclust"],
                                       cluster_score_method="ARI")

        elapsed = time.perf_counter() - start

        print("slice {} repeat {}: {:.2f} s".format(section_id, repeat + 1, elapsed))
        with open(os.path.join(dir_output,
                               "stcluster_runtime_{}.txt".format(section_id)),
                  "a+") as fp:
            fp.write(str(elapsed) + "\n")
