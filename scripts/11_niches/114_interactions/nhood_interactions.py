#%%
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import athena as ath



#%%
import sys
sys.path.append('/users/mensmeng/workspace/PCA_NHOODs_clean/workflow')
from utils.tools import nhood_filtering, prepare_info_for_filtering, leiden_clustering_scanpy
#################### DATA SETUP ####################


#%%
r=32
graph_type = "radius"


base_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/")
anndata_dir = base_dir / "anndatas"
cell_data_path = base_dir / "cell_metadata.parquet"
metadata_path = base_dir / "metadata.parquet"


#%%
# %%
path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood/Paper_CellCellNeighborhoods/low_res_merged/clusters_nhood_group_annotated.parquet')
orig_clusters = pd.read_parquet(path, engine='fastparquet')
print(f"Original clusters shape: {orig_clusters.shape}")



#%%

samples = [f for f in os.listdir(anndata_dir) if f.endswith('.pkl')]
## strip file extension
samples = [f[:-4] for f in samples]
# %%
sample_id = samples[55]
sample_path = anndata_dir / f"{sample_id}.pkl"
# %%
import pickle
with open(sample_path, 'rb') as f:
    adata = pickle.load(f)
print(adata)


# %%
## subset based on index
nhood_sample = orig_clusters.xs(sample_id, level='sample_id', drop_level=True)
nhood_sample = nhood_sample[['nhood', 'main_nhood']]
nhood_sample.index = nhood_sample.index.astype(str)
adata.obs = adata.obs.join(nhood_sample, how='inner')
print(adata)

#%%
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
import athena as ath

ath.pl.spatial(ad=adata, attr='label', edges=True, graph_key='radius_32', ax=axs[0], cbar=True)
ath.pl.spatial(ad=adata, attr='nhood', edges=True, graph_key='radius_32', ax=axs[1], cbar=True)
ath.pl.spatial(ad=adata, attr='main_nhood', edges=True, graph_key='radius_32', ax=axs[2], cbar=True)

# %%
nhood_type = 'tumorERG+_peritumoralCAF1(CD105+)'
adata_sub = adata[adata.obs['nhood'] == nhood_type]

#%%
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
ath.pl.spatial(ad=adata_sub, attr='nhood', edges=True, graph_key='radius_32', ax=axs[0], cbar=True)
ath.pl.spatial(ad=adata_sub, attr='label', edges=True, graph_key='radius_32', ax=axs[1], cbar=True)
# %%

save_dir = Path("/users/mensmeng/workspace/nhoods/PCa_NHood/interaction") / "nhood_annotated"
save_dir.mkdir(parents=True, exist_ok=True)
nhood_type = 'tumorERG+_peritumoralCAF1(CD105+)'

#%%
for nhood_type in orig_clusters['nhood'].unique():
    print(nhood_type)
    

    nhood_interaction_list = []
    from tqdm import tqdm
    for sample_id in tqdm(samples):
        sample_path = anndata_dir / f"{sample_id}.pkl"
        
        with open(sample_path, 'rb') as f:
            adata = pickle.load(f)
        #print(adata)


        ## subset based on index
        nhood_sample = orig_clusters.xs(sample_id, level='sample_id', drop_level=True)
        nhood_sample = nhood_sample[['nhood', 'main_nhood']]
        nhood_sample.index = nhood_sample.index.astype(str)
        adata.obs = adata.obs.join(nhood_sample, how='inner')
        #print(adata)


        ## check if nhood_type is in adata
        if nhood_type not in adata.obs['nhood'].values:
            print(f"{nhood_type} not in {sample_id}, skipping")
            continue
        if adata.obs['nhood'].value_counts()[nhood_type] < 20:
            print(f"{nhood_type} has less than 20 cells in {sample_id}, skipping")
            continue

        ## subset adata to only include cells in the specified neighborhood


        adata_sub = adata[adata.obs['nhood'] == nhood_type]

        ath.neigh.interactions(ad=adata_sub, attr='label', mode='proportion', prediction_type='observation',
                                graph_key='radius_32')
        ath.neigh.interactions(ad=adata, attr='label', mode='proportion', prediction_type='observation',
                                graph_key='radius_32')
        adata_exclude = adata[~adata.obs['nhood'].isin([nhood_type])]
        ath.neigh.interactions(ad=adata_exclude, attr='label', mode='proportion', prediction_type='observation',
                                graph_key='radius_32')
        


        key = "interaction_label_proportion_observation_radius_32"
        df_full = adata_exclude.uns[key]
        df_nhood = adata_sub.uns[key]

        ## rename columns
        df_full = df_full.rename(columns={'score': 'score_ROI'})
        df_nhood = df_nhood.rename(columns={'score': 'score_nhood'})

        df = df_full.join(df_nhood, how='outer')
        df = df.fillna(0)

        df['log2fc'] = np.log2((df['score_nhood'] + 1e-9) / (df['score_ROI'] + 1e-9))
        df['diff'] = df['score_nhood'] - df['score_ROI']
        ### extract log2fc and index and make a heatmap
        df.reset_index(inplace=True)
        nhood_interaction_list.append(df.assign(sample_id=sample_id))
        # df_heatmap = df.pivot(index='source_label', columns='target_label', values='log2fc')

        # plt.figure(figsize=(10, 8))
        # sns.heatmap(df_heatmap, cmap='coolwarm', center=0, annot=True, fmt=".2f", vmax=5, vmin=-5, annot_kws={"size": 6})
        # plt.title(f'Log2 Fold Change in Interaction Proportions\nfor {nhood_type} vs Full Sample')
        # plt.tight_layout()
        # plt.show()


    ######## plot the mean log2fc across all samples #######
    df_results = pd.concat(nhood_interaction_list)
    df_mean_log2fc = df_results.groupby(['source_label', 'target_label'])['log2fc'].mean().reset_index()
    df_heatmap = df_mean_log2fc.pivot(index='source_label', columns='target_label', values='log2fc')
    ### cap the values at -1 for better visualization
    df_heatmap = df_heatmap.clip(lower=-2, upper=np.inf)
    plt.figure(figsize=(15, 12))
    sns.heatmap(df_heatmap, cmap='coolwarm', center=0, annot=True, fmt=".2f", vmax=5, vmin=-5, annot_kws={"size": 4})
    plt.title(f'Log2 Fold Change in Interaction Proportions\nfor {nhood_type} vs Full Sample')
    plt.tight_layout()
    plot_dir = save_dir / "log2fc"
    plot_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(plot_dir / f"{nhood_type}_log2fc_heatmap.png", dpi=300)
    plt.show()

    ######## plot the mean diff across all samples #######
    df_results = pd.concat(nhood_interaction_list)
    df_mean_diff = df_results.groupby(['source_label', 'target_label'])['diff'].mean().reset_index()
    df_heatmap = df_mean_diff.pivot(index='source_label', columns='target_label', values='diff')
    ### cap the values at -1 for better visualization
    #df_heatmap = df_heatmap.clip(lower=-2, upper=np.inf)
    plt.figure(figsize=(15, 12))
    sns.heatmap(df_heatmap, cmap='coolwarm', center=0, annot=True, fmt=".2f", vmax=.5, vmin=-.5, annot_kws={"size": 4})
    plt.title(f'Difference in Interaction Proportions\nfor {nhood_type} vs Full Sample')
    plot_dir = save_dir / "diff"
    plot_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(plot_dir / f"{nhood_type}_diff_heatmap.png", dpi=300)
    plt.tight_layout()
    plt.show()



# %%
from matplotlib.colors import Normalize

norm = Normalize(0, 1)
ath.pl.interactions(ad=adata_sub, attr='label', mode='proportion',
                    prediction_type='observation', graph_key='radius_32',
                    ax=None, norm=norm, cmap='Blues', cbar=True)
### smaller fontsize
plt.xticks(fontsize=4)
plt.yticks(fontsize=4)
plt.title("focusing on " + nhood_type)
plt.tight_layout()
#

from matplotlib.colors import Normalize

norm = Normalize(0, 1)
ath.pl.interactions(ad=adata_exclude, attr='label', mode='proportion',
                    prediction_type='observation', graph_key='radius_32',
                    ax=None, norm=norm, cmap='Blues', cbar=True)
### smaller fontsize
plt.xticks(fontsize=4)
plt.yticks(fontsize=4)
plt.title("excluding " + nhood_type)
plt.tight_layout()

ath.pl.interactions(ad=adata, attr='label', mode='proportion',
                    prediction_type='observation', graph_key='radius_32',
                    ax=None, norm=norm, cmap='Blues', cbar=True)
### smaller fontsize
plt.xticks(fontsize=4)
plt.yticks(fontsize=4)
plt.title("full sample")
plt.tight_layout()

# %%
### 
adj = adata.obsp['radius_32'].toarray()
## check if adj is symmetric
print(np.allclose(adj, adj.T))
# %%


## 

# %%
