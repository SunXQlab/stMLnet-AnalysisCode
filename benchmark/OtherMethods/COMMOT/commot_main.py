import os
import time
import psutil
import csv
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import commot as ct

current_dir = os.getcwd()
print("current path:", current_dir)

def show_info(start):
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    memory = info.uss/1024/1024/1024
    return memory

start = show_info('strat')
start_time = time.time()

res_path = "E:/AA-luluyan-phd/code/01_cell_cell_inteaction/stMLnet/apply_in_stBC/COMMOT/result/"

# load st expression, spatial location and deconvolution annotation information
data_path = "E:/AA-luluyan-phd/code/01_cell_cell_inteaction/stMLnet/apply_in_stBC/COMMOT/input/"
st_count = pd.read_csv(data_path + "st_count.csv").to_numpy().astype(float).T
gene_id = pd.read_csv(data_path + "gene_ID.csv")
sp_anno = pd.read_csv(data_path + "spot_anno.csv").to_numpy()
spatial_locs = pd.read_csv(data_path + "spatial.locs.csv").to_numpy().astype(float)

gene = gene_id.iloc[:,0].tolist()
var = pd.DataFrame(index=gene)
adata = anndata.AnnData(X=st_count, var=var)
adata.obsm['spatial'] = spatial_locs
adata.obs['celltype'] = sp_anno

# Preprocessing the data
adata.var_names_make_unique()
adata.raw = adata
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

adata_dis500 = adata.copy()
print(adata_dis500)

# spatial communication inference
# use CellChatDB ligand-receptor database
df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
print('the prior cellchat database is',df_cellchat.shape)  # (1209, 4)

df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)
print('the filtered lr pair is',df_cellchat_filtered.shape)
print(df_cellchat_filtered.head())  # (300, 4)

ct.tl.spatial_communication(adata_dis500, database_name='cellchat',
                            df_ligrec=df_cellchat_filtered,
                            dis_thr=500,
                            heteromeric=True,
                            pathway_sum=True) 

ct.tl.communication_direction(adata_dis500, database_name='cellchat', pathway_name=None, k=5)
adata_dis500.write(res_path +'adata_pw.h5ad')
lr_keys = [s.split('-')[2:] for s in list(adata_dis500.obsp.keys())]

result = []
for lr in lr_keys:
    print('calculate the communication score of',lr)
    ct.tl.cluster_communication(adata_dis500, lr_pair=lr, database_name='cellchat', pathway_name=None,
                                clustering='celltype',
                                n_permutations=100)
    comm_mtx = adata_dis500.uns['commot_cluster-celltype-cellchat-'+lr[0]+'-'+lr[1]]['communication_matrix']
    comm_mtx = comm_mtx.reset_index()
    comm_mtx.rename(columns={'index': 'Sender'}, inplace=True)
    comm_mtx_df = comm_mtx.melt(id_vars='Sender', var_name='Receiver', value_name='score')
    comm_mtx_df['Ligand'] = lr[0]
    comm_mtx_df['Receptor'] = lr[1]
    result.append(comm_mtx_df)

result = pd.concat(result, ignore_index=True)
end = show_info('end')
end_time = time.time()
run_time = (end_time - start_time) / 60
print(f"Training time is: {run_time} mins")
print('total memory used '+str(end-start) + 'GB')
adata_dis500.write(res_path +'adata_pw_new.h5ad')
result.to_csv(res_path+'result.csv', index=False)

summary_name = 'commot-'+'cellchat'+'-sum-'+'receiver'
summary_abrv = 'r'
lr_pair: tuple = ('total','total')
comm_sum = adata_dis500.obsm[summary_name][summary_abrv+'-'+lr_pair[0]+'-'+lr_pair[1]].values.reshape(-1,1)
cell_weight = np.ones_like(comm_sum).reshape(-1,1)

np.savetxt(data_path+'/comm_sum_pw.csv', comm_sum, delimiter=',', fmt='%d')
np.savetxt(data_path+'/cell_weight_pw.csv', cell_weight, delimiter=',', fmt='%d')


