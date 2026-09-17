"""
stCluster on the mHypothalamus (MERFISH mouse hypothalamus) dataset.

Adapted from https://github.com/hannshu/stCluster

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download mHypothalamus data from https://zenodo.org/records/10698909

Expected layout:
  mHypothalamus/MERFISH_Animal1_cnts.xlsx   (one sheet per section)
  mHypothalamus/MERFISH_Animal1_info.xlsx   (same sheet names)

Outputs written to dir_output:
  stcluster_mhypo_<section>.csv             mclust labels (cell, mclust)
  stcluster_embeddings_mhypo_<section>.csv  stCluster embedding, cells x 30
  stcluster_aris_mhypo.txt                  one "mHypo<section> <ARI>" line per run
"""

import os
import warnings

import numpy as np
import scanpy as sc

from st_loading_utils import load_mHypothalamus
from stCluster.run import evaluate_embedding
from stCluster.train import train

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run stCluster for a selected section index (from 1 to 5)
section_index = 1
# ---------------------------------------------------------------------------

sections = ["-0.04", "-0.09", "-0.14", "-0.19", "-0.24"]

section_id = sections[section_index - 1]
n_cluster = 8

os.makedirs(dir_output, exist_ok=True)

ad = load_mHypothalamus(root_dir=dir_input, section_id=section_id)

# Preprocessing. No highly_variable_genes step because the dataset has less than 3000 genes
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
    os.path.join(dir_output, "stcluster_mhypo_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "stcluster_embeddings_mhypo_{}.csv".format(section_id)),
           np.asarray(ad.obsm["embedding"]), delimiter=",")
with open(os.path.join(dir_output, "stcluster_aris_mhypo.txt"), "a+") as fp:
    fp.write("mHypo" + section_id + " " + str(ARI) + "\n")
