"""
GraphST on the DLPFC12 dataset.

Adapted from https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the BenchmarkST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download DLPFC12 data from https://zenodo.org/records/10698880

Outputs written to dir_output:
  graphst_slice_<slice>.csv                final domain labels (barcode, domain)
  graphst_embeddings_slice_<slice>.csv     GraphST embedding
  graphst_embeddings_pca_slice_<slice>.csv 20-PC embedding mclust was run on
  graphst_aris.txt                         one "DLPFC<slice> <ARI>" line appended per run
"""

import os
import warnings

import numpy as np
import pandas as pd
import torch
from sklearn import metrics

from GraphST import GraphST
from GraphST.utils import clustering
from st_loading_utils import load_DLPFC

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run GraphST for a selected slice index (from 1 to 12)
slice_index = 9
# ---------------------------------------------------------------------------

device_name = "auto"

slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

section_id = str(slices[slice_index - 1])
n_clusters = 5 if 5 <= slice_index <= 8 else 7

radius = 50
tool = "mclust"

if device_name == "auto":
    device_name = "cuda" if torch.cuda.is_available() else "cpu"
device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

ad = load_DLPFC(root_dir=dir_input, section_id=section_id)

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
    os.path.join(dir_output, "graphst_slice_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output, "graphst_embeddings_slice_{}.csv".format(section_id)),
           np.asarray(ad.obsm["emb"]), delimiter=",")
np.savetxt(os.path.join(dir_output,
                        "graphst_embeddings_pca_slice_{}.csv".format(section_id)),
           np.asarray(ad.obsm["emb_pca"]), delimiter=",")
with open(os.path.join(dir_output, "graphst_aris.txt"), "a+") as fp:
    fp.write("DLPFC" + section_id + " " + str(ARI) + "\n")
