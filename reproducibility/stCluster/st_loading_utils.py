import os

import anndata
import numpy as np
import pandas as pd
import scanpy as sc


# for loading DLPFC12 data
def load_DLPFC(root_dir='../benchmarking_data/DLPFC12', section_id='151507'):
    # 151507, ..., 151676 12 in total
    ad = sc.read_visium(path=os.path.join(root_dir, section_id),
                        count_file=section_id + '_filtered_feature_bc_matrix.h5')
    ad.var_names_make_unique()

    gt_dir = os.path.join(root_dir, section_id, 'gt')
    gt_df = pd.read_csv(os.path.join(gt_dir, 'tissue_positions_list_GTs.txt'),
                        sep=',', header=None, index_col=0)
    ad.obs['cluster'] = gt_df.loc[:, 6]
    keep_bcs = ad.obs.dropna().index
    ad = ad[keep_bcs].copy()
    ad.obs['cluster'] = ad.obs['cluster'].astype(int).astype(str)
    return ad


# for loading BC data
# cluster = 20
def load_BC(root_dir='../benchmarking_data/BC', section_id='section1'):
    # section1
    ad = sc.read_visium(path=os.path.join(root_dir, section_id),
                        count_file=section_id + '_filtered_feature_bc_matrix.h5')
    ad.var_names_make_unique()

    gt_dir = os.path.join(root_dir, section_id, 'gt')
    gt_df = pd.read_csv(os.path.join(gt_dir, 'tissue_positions_list_GTs.txt'),
                        sep=',', header=None, index_col=0)
    ad.obs['cluster'] = gt_df.loc[:, 6].astype(int)
    ad.obs['cluster'] += 1
    keep_bcs = ad.obs.dropna().index
    ad = ad[keep_bcs].copy()
    ad.obs['cluster'] = ad.obs['cluster'].astype(int).astype(str)
    return ad


# for loading mMAMP (mouse brain anterior) data
# cluster = 52
def load_mMAMP(root_dir='../benchmarking_data/mMAMP', section_id='MA'):
    ad = sc.read_visium(path=os.path.join(root_dir, section_id),
                        count_file=section_id + '_filtered_feature_bc_matrix.h5')
    ad.var_names_make_unique()

    gt_dir = os.path.join(root_dir, section_id, 'gt')
    gt_df = pd.read_csv(os.path.join(gt_dir, 'tissue_positions_list_GTs.txt'),
                        sep='\t', header=0, index_col=0)
    ad.obs = gt_df
    ad.obs['cluster'] = ad.obs['ground_truth']
    return ad


# for loading mHypothalamus data
# already preprocessed? Xs are floats
def load_mHypothalamus(root_dir='../benchmarking_data/mHypothalamus', section_id='0.26'):
    # section id = '0.26', '0.21', '0.16', '0.11', '0.06', '0.01', '-0.04', '-0.09',
    #              '-0.14', '-0.19', '-0.24', '-0.29' 12 in total
    info_file = os.path.join(root_dir, 'MERFISH_Animal1_info.xlsx')
    cnts_file = os.path.join(root_dir, 'MERFISH_Animal1_cnts.xlsx')
    xls_cnts = pd.ExcelFile(cnts_file)
    df_cnts = pd.read_excel(xls_cnts, section_id)

    xls_info = pd.ExcelFile(info_file)
    df_info = pd.read_excel(xls_info, section_id)

    obs_ = df_info
    if len(df_info.columns) == 5:
        obs_.columns = ['psuedo_barcodes', 'x', 'y', 'cluster', 'Neuron_cluster_ID']
    elif len(df_info.columns) == 6:
        obs_.columns = ['psuedo_barcodes', 'x', 'y', 'cell_types',
                        'Neuron_cluster_ID', 'cluster']
    obs_.index = obs_['psuedo_barcodes'].tolist()

    var_ = df_cnts.iloc[:, 0]
    var_ = pd.DataFrame(var_)

    ad = anndata.AnnData(X=df_cnts.iloc[:, 1:].T, obs=obs_, var=var_)
    spatial = np.vstack((ad.obs['x'].to_numpy(), ad.obs['y'].to_numpy()))
    ad.obsm['spatial'] = spatial.T
    return ad
