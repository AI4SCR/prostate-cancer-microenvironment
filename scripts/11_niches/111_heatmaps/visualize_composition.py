#%%
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np
sys.path.append("/users/mensmeng/workspace/PCA_NHOODs_clean/robustness")


#%%
output_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/visualization/composition")
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Results will be saved to {output_dir}")

#%%
cluster_path = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/annotation/clusters_annotated_v2.parquet")
df_clusters = pd.read_parquet(cluster_path, engine='fastparquet')
print(f"Clusters shape: {df_clusters.shape}")
df_clusters.reset_index(inplace=True)






#%%
from utils.visualization import transform_for_heatmap, calculate_zscore, calculate_lfc, calculate_freqs, create_plot


nhood = 'niche'

df_clusters[nhood] = df_clusters[nhood].astype('category')
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

#show plot
plt.show()
#%%
plot_path = output_dir / f"{nhood}_heatmap.png"
p.savefig(plot_path, dpi=300, bbox_inches='tight')
heatmap_data.to_parquet(output_dir / f"{nhood}_heatmap_data.parquet", engine="fastparquet")
print(f"Saved heatmap to {plot_path}")


# %%
df_abundance = df_clusters[['sample_id', 'niche']].copy()
# %%
df_abundance['sample_id'] = df_abundance['sample_id'].astype('category')
df_abundance['niche'] = df_abundance['niche'].astype('category')

#%%
df_counts = df_abundance.groupby(['sample_id', 'niche']).size().reset_index(name='count')
sample_sizes = df_abundance['sample_id'].value_counts().reset_index()
sample_sizes.columns = ['sample_id', 'total_cells']
df_counts = df_counts.merge(sample_sizes, on='sample_id')
# %%
filter_by_presence = True
if filter_by_presence:
    min_cells = 15
    ## remove all rows where a niche has less than min_cells in a sample
    df_counts = df_counts[df_counts['count'] >= min_cells]
    print(f"Filtered df_counts shape: {df_counts.shape}")
    df_num_samples = df_counts.groupby('niche')['sample_id'].nunique().reset_index()
    df_num_samples.columns = ['niche', 'num_samples']
# %%
df_counts['frequency'] = df_counts['count'] / df_counts['total_cells']
df_means = df_counts.groupby('niche')['frequency'].mean().reset_index()
df_means.columns = ['niche', 'mean_frequency']
df_means = df_means.sort_values(by='mean_frequency', ascending=False)
df_stds = df_counts.groupby('niche')['frequency'].std().reset_index()
df_medians = df_counts.groupby('niche')['frequency'].median().reset_index()
df_medians.columns = ['niche', 'median_frequency']
df_stds.columns = ['niche', 'std_frequency']
df_stats = df_means.merge(df_stds, on='niche')
df_stats = df_stats.merge(df_medians, on='niche')

# %%
df_stats = df_stats.merge(df_num_samples, on='niche')
# %%
df_stats.to_parquet(output_dir / "niche_abundance_stats.parquet", engine="fastparquet")
print(f"Saved niche abundance stats to {output_dir / 'niche_abundance_stats.parquet'}")
# %%
