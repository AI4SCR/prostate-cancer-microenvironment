
#%%
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
import anndata as ad
import athena as ath
import pickle
sys.path.append("/users/mensmeng/workspace/PCA_NHOODs_clean/robustness")

#%% 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--sample', type=str, required=True, help='Sample ID to process')
args = parser.parse_args()
sample = args.sample
print(f"Processing sample: {sample}")



#%%
output_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/visualization/interactions")
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Results will be saved to {output_dir}")

#%%
cluster_path = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/annotation/clusters_annotated.parquet")
df_clusters = pd.read_parquet(cluster_path, engine='fastparquet')
print(f"Clusters shape: {df_clusters.shape}")



# %%
anndata_dir = Path('/users/mensmeng/workspace/nhoods/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas')

attr_list = ['label', 'niche', 'meta_niche']
sample_list = df_clusters.index.get_level_values('sample_id').unique().tolist()


############ compute interactions ##################
#%%
cell_types = df_clusters['label'].unique().tolist()
## sort alphabetically
cell_types.sort()
interaction_matrix = pd.DataFrame(0, index=cell_types, columns=cell_types)
graph_key = 'radius_32'

niche_ids = df_clusters['niche'].unique().tolist()
niche_ids.sort()

anndata_save_dir = output_dir / "anndatas_with_interactions"
anndata_save_dir.mkdir(parents=True, exist_ok=True)

#%%
from utils_interactions import compute_interaction_matrices, plot_interaction_matrix
#%%

print(f"Processing sample: {sample}")
anndata_path = anndata_dir / f"{sample}.pkl"
if not anndata_path.exists():
    print(f"Anndata file {anndata_path} does not exist")
with open(anndata_path, 'rb') as f:
    adata = pickle.load(f)
obs = adata.obs
sample_data = df_clusters.loc[sample]
print(f'shape of sample data: {sample_data.shape} ')
print(f'shape of obs: {obs.shape} ')
## make indices to string
sample_data.index = sample_data.index.astype(str)

#assert obs.index in data.index, f"Anndata index {obs.index} does not match data index {data.index}"
obs = obs.merge(sample_data, left_index=True, right_index=True, how='left', suffixes=('', '_full'))
#obs = obs.reset_index('sample_id')
adata.obs = obs

interaction_freqs = compute_interaction_matrices(adata, cell_types=cell_types, graph_key=graph_key)
title = f"Cell-Cell Interaction Matrix for Sample {sample}"
fig = plot_interaction_matrix(interaction_freqs, title=title)
fig_path = output_dir / sample / "interaction_overall.png"
fig_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(fig_path, dpi=300)

adata.uns['interaction_freqs_overall'] = interaction_freqs

for niche in niche_ids:
    print(f"Processing niche: {niche}")
    niche_adata = adata[adata.obs['niche'] == niche].copy()
    if niche_adata.n_obs <= 15:
        print(f"Skipping niche {niche} with only {niche_adata.n_obs} cells")
        continue
    interaction_freqs = compute_interaction_matrices(niche_adata, cell_types=cell_types, graph_key=graph_key)
    adata.uns[f'interaction_freqs_niche_{niche}'] = interaction_freqs
    title = f"Cell-Cell Interaction Matrix for Sample {sample}, Niche {niche}"
    fig = plot_interaction_matrix(interaction_freqs, title=title)
    fig_path = output_dir / sample / f"interaction_niche_{niche}.png"
    fig.savefig(fig_path, dpi=300)

# save updated anndata with interaction matrices
anndata_save_path = anndata_save_dir / f"{sample}_with_interactions.pkl"
anndata_save_path.parent.mkdir(parents=True, exist_ok=True)
with open(anndata_save_path, 'wb') as f:
    pickle.dump(adata, f)
    


#%%
# ## write sample_ids to txt
# sample_ids_path = Path("/users/mensmeng/workspace/nhoods") / "sample_ids.txt"
# with open(sample_ids_path, 'w') as f:   
#     for sample in sample_list:
#         f.write(f"{sample}\n")
# print(f"Wrote processed sample IDs to {sample_ids_path}")# %%
# # %%
