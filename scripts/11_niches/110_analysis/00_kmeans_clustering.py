#%%
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

import sys
sys.path.append("/users/mensmeng/workspace/nhoods/PCa/05_nhoods/PCA_NHOODs_clean/robustness")



###  DATA SETUP ####################
r=32
graph_type = "radius"
#data_type = "freq"  # "count" or "freq"
data_type = "freq"
assert data_type in ["count", "freq"], f"data_type must be 'count' or 'freq', got {data_type}"


count_base_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood/final_analysis/evaluation/count/CellCellNeighborhoods")
data_path = count_base_dir / f"graph_type={graph_type}-radius={r}" / "data.parquet"
cell_data_path = count_base_dir / "cell_metadata.parquet"
metadata_path = count_base_dir / "metadata.parquet"

data = pd.read_parquet(data_path, engine="fastparquet")
cell_metadata = pd.read_parquet(cell_data_path, engine="fastparquet")
metadata = pd.read_parquet(metadata_path, engine="fastparquet")


### DATA FILTERING #####################
threshold = 2

data = pd.read_parquet(data_path)
data_filt = data[data.sum(axis=1) > threshold]

print(f"Filtered out {data.shape[0] - data_filt.shape[0]} cells with less than {threshold+1} neighbors")
data_orig = data.copy()
data = data_filt.copy()

######### FREQUENCY VS. COUNT ##########
if data_type == "count":
    print("Using count data for clustering")
elif data_type == "freq":
    print("Using frequency data for clustering")
    data = data.div(data.sum(axis=1), axis=0)  # normalize to percentage




### CLUSTERING PREP ###
#%%
df_results = pd.DataFrame(index=data.index)
k = 24
seed = 686

output_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/clustering")
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Results will be saved to {output_dir}")


#%%
from utils.clustering import perform_kmeans_clustering
print(f"Running k-means with k={k}, seed={seed}")
cluster_data, cluster_name = perform_kmeans_clustering(data, k=k, random_state=seed)


print(f"Cluster data shape: {cluster_data.shape}")
print(f"Cluster name: {cluster_name}")

from utils.clustering import wrapper_nhood_filtering
cluster_data = wrapper_nhood_filtering(cluster_data, cluster_name, metadata, cell_metadata, min_cells=15, min_patients=5)
print(f"Cluster data shape after filtering: {cluster_data.shape}")  

final_cluster_name = cluster_data.columns[-1]
df_results = df_results.join(cluster_data[[final_cluster_name]], how='inner')
assert data.index.equals(df_results.index), f"data index {data.index} does not match df_results index {df_results.index}"
assert data.shape[0] == df_results.shape[0], f"data shape {data.shape} does not match df_results shape {df_results.shape}"
print(f"df_results shape: {df_results.shape}")


###########VISUALIZATION ##############
# %%


from utils.visualization import transform_for_heatmap, calculate_zscore, calculate_lfc, calculate_freqs, create_plot

df_clusters = df_results.join(cell_metadata, how='right')
df_clusters.fillna('unassigned', inplace=True)
# assert df_clusters.shape[0] == df_results.shape[0], f"df_clusters shape {df_clusters.shape} does not match df_results shape {df_results.shape}"
df_clusters.reset_index(inplace=True)
nhood = final_cluster_name



df_clusters[final_cluster_name] = df_clusters[nhood].astype('category')
df_clusters['label'] = df_clusters['label'].astype('category')
df_counts, df_freqs =  calculate_freqs(df_clusters, cell_type_col='label', nhood_col=nhood, sample_col='sample_id')

target='z_score_sample_normalized'
metric='Frequency normalized sample-wise (z-score)'
x='mean_nhood_celltype_sample'
m='mean_celltype_sample'
s='std_celltype_sample'

freq_df = calculate_zscore(df=df_freqs,
x_col=x, mean_col=m, std_col=s, colname=target)



heatmap_data = transform_for_heatmap(freq_df, target_col=target,
                                cell_type_col='label', nhood_col=nhood)

p = create_plot(df=heatmap_data, 
                    metric=metric,
                    nhood=nhood,
                    color_scheme='coolwarm',
                    upper_limit=3, lower_limit=-3,
                    cluster_celltypes=False, cluster_neighborhoods=True)
print(f"Saving heatmap to {output_dir / f'{nhood}_heatmap.png'}")
p.savefig(output_dir / f"{nhood}_heatmap.png", dpi=300, bbox_inches='tight')


#%%
df_clusters.set_index(['sample_id', 'object_id'], inplace=True)
#%%
df_clusters.to_parquet(output_dir / "clusters.parquet", engine="fastparquet")
print(f"Saved clustering results to {output_dir / 'clusters.parquet'}")

# %%
