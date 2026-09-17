"""
SpaGCN on the mHypothalamus (MERFISH mouse hypothalamus) dataset.

Adapted from https://github.com/jianhuupenn/SpaGCN
and https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download mHypothalamus data from https://zenodo.org/records/10698909

Expected layout:
  mHypothalamus/MERFISH_Animal1_cnts.xlsx   (one sheet per section)
  mHypothalamus/MERFISH_Animal1_info.xlsx   (same sheet names)

Outputs written to dir_output:
  spagcn_mHypothalamus_<section>.csv              refined labels (cell, refined_pred)
  spagcn_embeddings_mHypothalamus_<section>.csv   SpaGCN latent embedding, cells x 50
  spagcn_aris_mHypothalamus.txt                   one "mHypothalamus<section> <ARI>" line
"""

import os
import random
import warnings

import numpy as np
import scanpy as sc
import torch
from sklearn.metrics import adjusted_rand_score

import SpaGCN as spg
from st_loading_utils import load_mHypothalamus

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run SpaGCN for a selected section index (from 1 to 5)
section_index = 1
# ---------------------------------------------------------------------------

sections = ["-0.04", "-0.09", "-0.14", "-0.19", "-0.24"]

section_id = sections[section_index - 1]
n_clusters = 8

r_seed = t_seed = n_seed = 100

os.makedirs(dir_output, exist_ok=True)

adata = load_mHypothalamus(root_dir=dir_input, section_id=section_id)
x_array = adata.obs["x"].tolist()
y_array = adata.obs["y"].tolist()
x_pixel = x_array
y_pixel = y_array

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
    os.path.join(dir_output, "spagcn_mHypothalamus_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "spagcn_embeddings_mHypothalamus_{}.csv".format(section_id)),
           np.asarray(adata.obsm["X_spagcn"]), delimiter=",")
with open(os.path.join(dir_output, "spagcn_aris_mHypothalamus.txt"), "a+") as fp:
    fp.write("mHypothalamus" + section_id + " " + str(ARI) + "\n")
