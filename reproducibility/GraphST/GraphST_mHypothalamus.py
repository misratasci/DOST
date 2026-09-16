"""
GraphST on the mHypothalamus (MERFISH mouse hypothalamus) dataset.

Adapted from https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download mHypothalamus data from https://zenodo.org/records/10698909

Expected layout:
  mHypothalamus/MERFISH_Animal1_cnts.xlsx   (one sheet per section)
  mHypothalamus/MERFISH_Animal1_info.xlsx   (same sheet names)

Outputs written to dir_output:
  graphst_mHypothalamus_<section>.csv                 domain labels (cell, domain)
  graphst_embeddings_mHypothalamus_<section>.csv      GraphST embedding
  graphst_embeddings_pca_mHypothalamus_<section>.csv  20-PC embedding mclust used
  graphst_aris_mHypothalamus.txt                      one "mHypothalamus<section> <ARI>" line
"""

import os
import warnings

import numpy as np
import pandas as pd
import torch
from sklearn import metrics

from GraphST import GraphST
from GraphST.utils import clustering
from st_loading_utils import load_mHypothalamus

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run GraphST for a selected section index (from 1 to 5)
section_index = 1
# ---------------------------------------------------------------------------

device_name = "auto"

# Section IDs (sheet names in the two xlsx files)
sections = ["-0.04", "-0.09", "-0.14", "-0.19", "-0.24"]

section_id = sections[section_index - 1]
n_clusters = 8

radius = 50
tool = "mclust"

if device_name == "auto":
    device_name = "cuda" if torch.cuda.is_available() else "cpu"
device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

ad = load_mHypothalamus(root_dir=dir_input, section_id=section_id)

model = GraphST.GraphST(ad, device=device)
ad = model.train()

if tool == "mclust":
    clustering(ad, n_clusters, radius=radius, method=tool, refinement=True)
elif tool in ["leiden", "louvain"]:
    clustering(ad, n_clusters, radius=radius, method=tool,
               start=0.1, end=2.0, increment=0.01, refinement=False)

ad = ad[~pd.isnull(ad.obs["original_clusters"])]

ARI = metrics.adjusted_rand_score(ad.obs["domain"], ad.obs["original_clusters"])
ad.uns["ARI"] = ARI

print("Dataset:", section_id)
print("ARI:", ARI)

ad.obs["domain"].to_csv(
    os.path.join(dir_output, "graphst_mHypothalamus_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "graphst_embeddings_mHypothalamus_{}.csv".format(section_id)),
           np.asarray(ad.obsm["emb"]), delimiter=",")
np.savetxt(os.path.join(dir_output,
                        "graphst_embeddings_pca_mHypothalamus_{}.csv".format(section_id)),
           np.asarray(ad.obsm["emb_pca"]), delimiter=",")
with open(os.path.join(dir_output, "graphst_aris_mHypothalamus.txt"), "a+") as fp:
    fp.write("mHypothalamus" + section_id + " " + str(ARI) + "\n")
