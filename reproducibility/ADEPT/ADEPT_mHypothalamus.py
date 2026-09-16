"""
ADEPT on the mHypothalamus (MERFISH mouse hypothalamus) dataset.

Adapted from https://benchmarkst-reproducibility.readthedocs.io/

We used the data links provided by the BenchmarkST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download mHypothalamus data from https://zenodo.org/records/10698909

Expected layout:
  mHypothalamus/MERFISH_Animal1_cnts.xlsx   (one sheet per section)
  mHypothalamus/MERFISH_Animal1_info.xlsx   (same sheet names)

Outputs written to dir.output:
  mHypothalamus_<section>_adept.txt              final cluster labels (one per cell)
  adept_embeddings_mHypothalamus_<section>.csv   final latent embedding
  adept_aris_mHypothalamus.txt                   one "mHypothalamus<section> <ARI>" line
"""

import os
import warnings
from argparse import Namespace

import numpy as np

import GAAE
from GAAE.utils import DE_num_calc, filter_num_calc, impute, initialize

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/mHypothalamus"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Run ADEPT for a selected section index (from 1 to 5)
section_index = 1
# ---------------------------------------------------------------------------

# Section IDs (sheet names in the two xlsx files)
sections = ["-0.04", "-0.09", "-0.14", "-0.19", "-0.24"]

section_id = sections[section_index - 1]
cluster_num = 8

# Hyperparameters set here are from
# https://benchmarkst-reproducibility.readthedocs.io/en/latest/ADEPT.html
args = Namespace(
    data_dir=dir_input,
    gt_dir=dir_input,
    input_data=section_id,
    cluster_num=cluster_num,
    impute_cluster_num=[cluster_num],
    radius=150,
    de_candidates="None",
    no_de=0,
    use_mean=0,
    impute_runs=2,
    runs=1,
    gt=1,
    use_hvgs=0,
    use_preprocessing=1,
    save_fig=0,
    filter_nzr=0.15,
    filter_num=None,
    de_nzr_min=0.299,
    de_nzr_max=0.399,
    use_gpu_id="0",
)

cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
cache_file = os.path.join(cache_dir, "mHypo" + section_id + ".txt")
os.makedirs(cache_dir, exist_ok=True)
os.makedirs(dir_output, exist_ok=True)

filter_num = filter_num_calc(args, args.filter_num)
print("optimized filter number = ", filter_num)

adata, adata_ori = initialize(args, filter_num)

if os.path.exists(cache_file):
    with open(cache_file, "r") as fp:
        de_top_k_list = [int(e) for e in fp.readlines()[0].strip().split(",")]
    print("previously cached de list = ", de_top_k_list)
else:
    de_top_k_list = DE_num_calc(args, adata)
    print("optimized de list = ", de_top_k_list)
    with open(cache_file, "a+") as fp:
        fp.write(",".join([str(i) for i in de_top_k_list]))

de_list_epoch = []
adata_list = []
if de_top_k_list != []:
    print("performing DEGs selection")
    for de_ in de_top_k_list:
        for cluster_n in args.impute_cluster_num:
            print("cluster_n = ", cluster_n)
            GAAE.get_kNN(adata, rad_cutoff=args.radius)
            ari_ini, ari_final, de_list, adata_out = GAAE.train_ADEPT_use_DE(
                adata, n_epochs=1000, num_cluster=int(cluster_n),
                dif_k=de_, device_id=args.use_gpu_id)
            de_list_epoch.append(de_list)
            adata_list.append(adata_out)
    g_union = set.union(*de_list_epoch)
    imputed_ad = impute(args, adata_list, g_union, de_top_k_list)
else:
    print("skip performing DEGs selection")
    imputed_ad = adata

GAAE.get_kNN(imputed_ad, rad_cutoff=args.radius)
ari_ini, ARI, de_list, adata_out = GAAE.train_ADEPT_use_DE(
    imputed_ad, n_epochs=1000, num_cluster=args.cluster_num,
    device_id=args.use_gpu_id)

print("Dataset:", section_id)
print("ARI:", ARI)

np.savetxt(os.path.join(dir_output, "mHypothalamus_{}_adept.txt".format(section_id)),
           adata_out.obs["mclust_impute"], fmt="%s")
np.savetxt(os.path.join(dir_output,
                        "adept_embeddings_mHypothalamus_{}.csv".format(section_id)),
           np.asarray(adata_out.obsm["ade_impute"]), delimiter=",")
with open(os.path.join(dir_output, "adept_aris_mHypothalamus.txt"), "a+") as fp:
    fp.write("mHypothalamus" + section_id + " " + str(ARI) + "\n")
