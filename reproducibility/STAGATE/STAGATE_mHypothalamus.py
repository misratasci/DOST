"""
STAGATE on the mHypothalamus (MERFISH mouse hypothalamus) dataset.

Adapted from https://github.com/QIFEIDKN/STAGATE_pyG
and https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download mHypothalamus data from https://zenodo.org/records/10698909

Expected layout:
  mHypothalamus/MERFISH_Animal1_cnts.xlsx   (one sheet per section)
  mHypothalamus/MERFISH_Animal1_info.xlsx   (same sheet names)

Outputs written to dir_output:
  stagate_mHypothalamus_<section>.csv             mclust labels (cell, mclust)
  stagate_embeddings_mHypothalamus_<section>.csv  STAGATE embedding, cells x 30
  stagate_aris_mHypothalamus.txt                  one "mHypothalamus<section> <ARI>" line
"""

import os
import warnings

import numpy as np
import scanpy as sc
import torch
from sklearn.metrics.cluster import adjusted_rand_score

import STAGATE_pyG as STAGATE
from st_loading_utils import load_mHypothalamus

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run STAGATE for a selected section index (from 1 to 5)
section_index = 1
# ---------------------------------------------------------------------------

device_name = "auto"

sections = ["-0.04", "-0.09", "-0.14", "-0.19", "-0.24"]

section_id = sections[section_index - 1]
n_clusters = 8

# STAGATE article recommends rad_cutoff chosen empirically
# so that each spot contains 6–15 neighbors on average
# for non 10x Visium data
rad_cutoff = 50 # 13.6 neighbors per cell on average for section_index=1

if device_name == "auto":
    device_name = "cuda:0" if torch.cuda.is_available() else "cpu"
device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

ad = load_mHypothalamus(root_dir=dir_input, section_id=section_id)

# No HVG selection because the MERFISH panel detects less than 3000 genes
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)

STAGATE.Cal_Spatial_Net(ad, rad_cutoff=rad_cutoff)
ad = STAGATE.train_STAGATE(ad, device=device)

sc.pp.neighbors(ad, use_rep="STAGATE")
sc.tl.umap(ad)
ad = STAGATE.mclust_R(ad, used_obsm="STAGATE", num_cluster=n_clusters)

obs_df = ad.obs.dropna()
ARI = adjusted_rand_score(obs_df["mclust"], obs_df["original_clusters"])

print("Dataset:", section_id)
print("ARI:", ARI)

ad.obs["mclust"].to_csv(
    os.path.join(dir_output, "stagate_mHypothalamus_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "stagate_embeddings_mHypothalamus_{}.csv".format(section_id)),
           np.asarray(ad.obsm["STAGATE"]), delimiter=",")
with open(os.path.join(dir_output, "stagate_aris_mHypothalamus.txt"), "a+") as fp:
    fp.write("mHypothalamus" + section_id + " " + str(ARI) + "\n")
