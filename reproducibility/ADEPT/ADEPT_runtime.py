"""
Runtime measurement for ADEPT on the DLPFC dataset.

Adapted from https://benchmarkst-reproducibility.readthedocs.io/

Times the full ADEPT pipeline repeated `n_repeats` times per slice. Nothing is saved except the timings.

We used the data links provided by the Benchmark ST study:
https://benchmarkst-reproducibility.readthedocs.io/en/latest/Data%20availability.html
Download DLPFC12 data from https://zenodo.org/records/10698880

Outputs written to dir.output:
  adept_runtime_<slice>.txt   one wall-clock time in seconds per line, appended
"""

import os
import time
import warnings
from argparse import Namespace

import GAAE
from GAAE.utils import DE_num_calc, filter_num_calc, impute, initialize

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Change this path to where you downloaded and unzipped the data
dir_input = "path/to/data/DLPFC12"

# Change this path to where you want to save the results
dir_output = "path/to/output/"

# Which slices to time (indices from 1 to 12) and how many repeats each
slice_indices = range(1, 13)
n_repeats = 5
# ---------------------------------------------------------------------------

# Slice IDs
slices = [151507, 151508, 151509, 151510,
          151669, 151670, 151671, 151672,
          151673, 151674, 151675, 151676]

cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache")
os.makedirs(cache_dir, exist_ok=True)
os.makedirs(dir_output, exist_ok=True)

for slice_index in slice_indices:
    section_id = str(slices[slice_index - 1])
    cluster_num = 5 if 5 <= slice_index <= 8 else 7

    # Same hyperparameters as ADEPT_DLPFC.py
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

    cache_file = os.path.join(cache_dir, "DLPFC" + section_id + ".txt")

    for repeat in range(n_repeats):
        start = time.perf_counter()

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

        GAAE.get_kNN(imputed_ad, rad_cutoff=args.radius)
        ari_ini, ARI, de_list, adata_out = GAAE.train_ADEPT_use_DE(
            imputed_ad, n_epochs=1000, num_cluster=args.cluster_num,
            device_id=args.use_gpu_id)

        elapsed = time.perf_counter() - start
        print("DLPFC{} repeat {}: {:.1f} s (ARI {:.4f})".format(
            section_id, repeat + 1, elapsed, ARI))

        with open(os.path.join(dir_output,
                               "adept_runtime_{}.txt".format(section_id)), "a+") as fp:
            fp.write(str(elapsed))
            fp.write("\n")
