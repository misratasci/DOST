"""
stCluster on the DLPFC12 dataset.

Adapted from https://github.com/hannshu/stCluster

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download DLPFC12 data from https://zenodo.org/records/10698880

Outputs written to dir_output:
  stcluster_slice_<slice>.csv             mclust labels (barcode, mclust)
  stcluster_embeddings_slice_<slice>.csv  stCluster embedding, spots x 30
  stcluster_aris.txt                      one "DLPFC<slice> <ARI>" line appended per run
"""

import os
import warnings

import numpy as np
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

# Run stCluster for a selected slice index (from 1 to 12)
slice_index = 9
# ---------------------------------------------------------------------------

slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

section_id = str(slices[slice_index - 1])
n_cluster = 5 if 5 <= slice_index <= 8 else 7

os.makedirs(dir_output, exist_ok=True)

ad = load_DLPFC(root_dir=dir_input, section_id=section_id)

# Preprocessing identical to stCluster's get_data function
sc.pp.highly_variable_genes(ad, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)

# Hyperparameters from https://stcluster.readthedocs.io/en/latest/section1.html
ad, g = train(ad, radius=150, ae_rate=0.8, adj_rate=0.2, pred_rate=0.3, seed=0)
print("embedding shape:", ad.obsm["embedding"].shape)

ad, score = evaluate_embedding(adata=ad, n_cluster=n_cluster,
                               cluster_method=["mclust"], cluster_score_method="ARI")
ARI = score["mclust"]

print("Dataset:", section_id)
print("ARI:", ARI)

ad.obs["mclust"].to_csv(
    os.path.join(dir_output, "stcluster_slice_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "stcluster_embeddings_slice_{}.csv".format(section_id)),
           np.asarray(ad.obsm["embedding"]), delimiter=",")
with open(os.path.join(dir_output, "stcluster_aris.txt"), "a+") as fp:
    fp.write("DLPFC" + section_id + " " + str(ARI) + "\n")
