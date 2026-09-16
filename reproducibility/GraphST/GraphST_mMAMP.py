"""
GraphST on the mMAMP (mouse brain anterior, "MA") dataset.

Adapted from https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download mMAMP data from https://zenodo.org/records/10698931

Expected layout:
  mMAMP/MA/MA_filtered_feature_bc_matrix.h5
  mMAMP/MA/spatial/
  mMAMP/MA/gt/tissue_positions_list_GTs.txt   (tab separated, header, 'ground_truth')

Outputs written to dir_output:
  graphst_mMAMP.csv                 final domain labels (barcode, domain)
  graphst_embeddings_mMAMP.csv      GraphST embedding (spots x 64)
  graphst_embeddings_pca_mMAMP.csv  20-PC embedding mclust was run on
  graphst_mMAMP_aris.txt            one "mABMA <ARI>" line appended per run
"""

import os
import warnings

import numpy as np
import pandas as pd
import torch
from sklearn import metrics

from GraphST import GraphST
from GraphST.utils import clustering
from st_loading_utils import load_mMAMP

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/mMAMP"

# Change this path to where you want to save the results
dir_output = "path/to/output/"
# ---------------------------------------------------------------------------

device_name = "auto"

section_id = "MA"
n_clusters = 52

radius = 50
tool = "mclust"

if device_name == "auto":
    device_name = "cuda" if torch.cuda.is_available() else "cpu"
device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

ad = load_mMAMP(root_dir=dir_input, section_id=section_id)

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

ad.obs["domain"].to_csv(os.path.join(dir_output, "graphst_mMAMP.csv"))
np.savetxt(os.path.join(dir_output, "graphst_embeddings_mMAMP.csv"),
           np.asarray(ad.obsm["emb"]), delimiter=",")
np.savetxt(os.path.join(dir_output, "graphst_embeddings_pca_mMAMP.csv"),
           np.asarray(ad.obsm["emb_pca"]), delimiter=",")
with open(os.path.join(dir_output, "graphst_mMAMP_aris.txt"), "a+") as fp:
    fp.write("mAB" + section_id + " " + str(ARI) + "\n")
