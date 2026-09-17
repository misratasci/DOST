"""
STAGATE on the DLPFC12 dataset.

Adapted from https://github.com/QIFEIDKN/STAGATE_pyG
and https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the BenchmarkST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download DLPFC12 data from https://zenodo.org/records/10698880

Outputs written to dir_output:
  stagate_slice_<slice>.csv             mclust labels (barcode, mclust)
  stagate_embeddings_slice_<slice>.csv  STAGATE embedding
  stagate_aris.txt                      one "DLPFC<slice> <ARI>" line appended per run
"""

import os
import warnings

import numpy as np
import scanpy as sc
import torch
from sklearn.metrics.cluster import adjusted_rand_score

import STAGATE_pyG as STAGATE
from st_loading_utils import load_DLPFC

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run STAGATE for a selected slice index (from 1 to 12)
slice_index = 9
# ---------------------------------------------------------------------------

device_name = "auto"

slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

section_id = str(slices[slice_index - 1])
n_clusters = 5 if 5 <= slice_index <= 8 else 7

rad_cutoff = 150

if device_name == "auto":
    device_name = "cuda:0" if torch.cuda.is_available() else "cpu"
device = torch.device(device_name)

os.makedirs(dir_output, exist_ok=True)

ad = load_DLPFC(root_dir=dir_input, section_id=section_id)
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

ad.obs["mclust"].to_csv(
    os.path.join(dir_output, "stagate_slice_{}.csv".format(section_id)))
np.savetxt(os.path.join(dir_output,
                        "stagate_embeddings_slice_{}.csv".format(section_id)),
           np.asarray(ad.obsm["STAGATE"]), delimiter=",")
with open(os.path.join(dir_output, "stagate_aris.txt"), "a+") as fp:
    fp.write("DLPFC" + section_id + " " + str(ARI) + "\n")
