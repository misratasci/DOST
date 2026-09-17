"""
Runtime measurement for SpaGCN on DLPFC12.

Times the full pipeline per slice: adjacency, gene filtering, normalization,
the search for l, the resolution search (which itself fits a model per candidate
resolution), training, prediction and label refinement.
Data loading and the coordinate extraction sit outside the timer.

Adapted from https://github.com/jianhuupenn/SpaGCN
and https://benchmarkst-reproducibility.readthedocs.io/

Outputs written to dir_output:
  spagcn_runtime_<slice>.txt   one wall-clock time in seconds per line, appended
"""

import os
import random
import time
import warnings

import numpy as np
import scanpy as sc
import torch

import SpaGCN as spg
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

slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

r_seed = t_seed = n_seed = 100

os.makedirs(dir_output, exist_ok=True)

for slice_index in slice_indices:
    section_id = str(slices[slice_index - 1])
    n_clusters = 5 if 5 <= slice_index <= 8 else 7

    for repeat in range(n_repeats):
        adata = load_DLPFC(root_dir=dir_input, section_id=section_id)
        x_array = adata.obs["array_row"].tolist()
        y_array = adata.obs["array_col"].tolist()
        x_pixel = adata.obsm["spatial"][:, 0].tolist()
        y_pixel = adata.obsm["spatial"][:, 1].tolist()

        start = time.perf_counter()

        adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel,
                                       x_pixel=x_pixel, y_pixel=y_pixel,
                                       image=None, beta=49, alpha=1, histology=False)
        spg.prefilter_genes(adata, min_cells=3)
        spg.prefilter_specialgenes(adata)
        sc.pp.normalize_per_cell(adata)
        sc.pp.log1p(adata)
        p = 0.5
        l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
        res = spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3,
                             lr=0.05, max_epochs=20,
                             r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
        clf = spg.SpaGCN()
        clf.set_l(l)
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)
        clf.train(adata, adj, init_spa=True, init="louvain", res=res,
                  tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob, adata = clf.predict(adata)
        adata.obs["pred"] = y_pred
        adata.obs["pred"] = adata.obs["pred"].astype("category")
        adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
        refined_pred = spg.refine(sample_id=adata.obs.index.tolist(),
                                  pred=adata.obs["pred"].tolist(),
                                  dis=adj_2d, shape="hexagon")
        adata.obs["refined_pred"] = refined_pred
        adata.obs["refined_pred"] = adata.obs["refined_pred"].astype("category")

        elapsed = time.perf_counter() - start

        print("slice {} repeat {}: {:.2f} s".format(section_id, repeat + 1, elapsed))
        with open(os.path.join(dir_output,
                               "spagcn_runtime_{}.txt".format(section_id)),
                  "a+") as fp:
            fp.write(str(elapsed) + "\n")
