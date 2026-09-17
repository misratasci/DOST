"""
STAGATE on the BC (human breast cancer, Visium) dataset.

Adapted from https://github.com/QIFEIDKN/STAGATE_pyG
and https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download BC data from https://zenodo.org/records/10698903

Expected layout:
  BC/section1/section1_filtered_feature_bc_matrix.h5
  BC/section1/spatial/
  BC/section1/gt/tissue_positions_list_GTs.txt

Outputs written to dir_output:
  stagate_bc.csv               mclust labels (barcode, mclust)
  stagate_embeddings_bc.csv    STAGATE embedding, spots x 30
  stagate_aris_bc.txt          one "HBRC1 <ARI>" line appended per run
"""

import os
import warnings

import numpy as np
import scanpy as sc
import torch
from sklearn.metrics.cluster import adjusted_rand_score

import STAGATE_pyG as STAGATE
from st_loading_utils import load_BC

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/BC"

# Change this path to where you want to save the results
dir_output = "path/to/output/"
# ---------------------------------------------------------------------------

device_name = "auto"

section_id = "section1"
n_clusters = 20
rad_cutoff = 450

if device_name == "auto":
    device_name = "cuda:0" if torch.cuda.is_available() else "cpu"
device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

ad = load_BC(root_dir=dir_input, section_id=section_id)
sc.pp.highly_variable_genes(ad, flavor="seurat_v3", n_top_genes=3000)
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

ad.obs["mclust"].to_csv(os.path.join(dir_output, "stagate_bc.csv"))
np.savetxt(os.path.join(dir_output, "stagate_embeddings_bc.csv"),
           np.asarray(ad.obsm["STAGATE"]), delimiter=",")
with open(os.path.join(dir_output, "stagate_aris_bc.txt"), "a+") as fp:
    fp.write("HBRC1 " + str(ARI) + "\n")
