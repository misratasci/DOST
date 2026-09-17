"""
SpaGCN on the DLPFC12 dataset.

Adapted from https://github.com/jianhuupenn/SpaGCN
and https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download DLPFC12 data from https://zenodo.org/records/10698880

Outputs written to dir_output:
  spagcn_slice_<slice>.csv              refined cluster labels (barcode, refined_pred)
  spagcn_embeddings_slice_<slice>.csv   SpaGCN latent embedding, spots x 50
  spagcn_aris.txt                       one "DLPFC<slice> <ARI>" line appended per run
"""

import os
import random
import warnings

import numpy as np
import scanpy as sc
import torch
from sklearn.metrics import adjusted_rand_score

import SpaGCN as spg
from st_loading_utils import load_DLPFC

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run SpaGCN for a selected slice index (from 1 to 12)
slice_index = 9
# ---------------------------------------------------------------------------

slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

section_id = str(slices[slice_index - 1])
n_clusters = 5 if 5 <= slice_index <= 8 else 7

r_seed = t_seed = n_seed = 100

os.makedirs(dir_output, exist_ok=True)

# 1. Load
adata = load_DLPFC(root_dir=dir_input, section_id=section_id)

x_array = adata.obs["array_row"].tolist()
y_array = adata.obs["array_col"].tolist()
x_pixel = adata.obsm["spatial"][:, 0].tolist()
y_pixel = adata.obsm["spatial"][:, 1].tolist()

# We set histology=False for fair comparison with other methods that do not use image info
adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel,
                               image=None, beta=49, alpha=1, histology=False)

spg.prefilter_genes(adata, min_cells=3)
spg.prefilter_specialgenes(adata)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

p = 0.5
l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
res = spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05,
                     max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

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

ARI = adjusted_rand_score(adata.obs["refined_pred"], adata.obs["original_clusters"])

print("Dataset:", section_id)
print("ARI:", ARI)

adata.obs["refined_pred"].to_csv(
    os.path.join(dir_output, "spagcn_slice_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "spagcn_embeddings_slice_{}.csv".format(section_id)),
           np.asarray(adata.obsm["X_spagcn"]), delimiter=",")
with open(os.path.join(dir_output, "spagcn_aris.txt"), "a+") as fp:
    fp.write("DLPFC" + section_id + " " + str(ARI) + "\n")
