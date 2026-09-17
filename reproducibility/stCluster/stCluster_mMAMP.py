"""
stCluster on the mMAMP (mouse brain anterior, "MA") dataset.

Adapted from https://github.com/hannshu/stCluster

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download mMAMP data from https://zenodo.org/records/10698931

Expected layout:
  mMAMP/MA/MA_filtered_feature_bc_matrix.h5
  mMAMP/MA/spatial/
  mMAMP/MA/gt/tissue_positions_list_GTs.txt   (tab separated, header, 'ground_truth')

Outputs written to dir_output:
  stcluster_mmamp_MA.csv             mclust labels (barcode, mclust)
  stcluster_embeddings_mmamp_MA.csv  stCluster embedding, spots x 30
  stcluster_aris_mmamp.txt           one "MA <ARI>" line appended per run
"""

import os
import warnings

import numpy as np
import scanpy as sc

from st_loading_utils import load_mMAMP
from stCluster.run import evaluate_embedding
from stCluster.train import train

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/mMAMP"

# Change this path to where you want to save the results
dir_output = "path/to/output/"
# ---------------------------------------------------------------------------

section_id = "MA"
n_cluster = 52

os.makedirs(dir_output, exist_ok=True)

ad = load_mMAMP(root_dir=dir_input, section_id=section_id)

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
    os.path.join(dir_output, "stcluster_mmamp_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "stcluster_embeddings_mmamp_{}.csv".format(section_id)),
           np.asarray(ad.obsm["embedding"]), delimiter=",")
with open(os.path.join(dir_output, "stcluster_aris_mmamp.txt"), "a+") as fp:
    fp.write(section_id + " " + str(ARI) + "\n")
