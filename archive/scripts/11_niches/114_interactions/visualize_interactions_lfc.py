
#%%
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import anndata as ad
import athena as ath
import pickle
sys.path.append("/users/mensmeng/workspace/PCA_NHOODs_clean/robustness")



#%%
output_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/visualization/interactions")
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Results will be saved to {output_dir}")

#%%
cluster_path = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/annotation/clusters_annotated.parquet")
df_clusters = pd.read_parquet(cluster_path, engine='fastparquet')
print(f"Clusters shape: {df_clusters.shape}")



# %%

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

anndata_dir = output_dir / "anndatas_with_interactions"

#%%
overall_interactions = []
overall_key = "interaction_freqs_overall"

for sample_id in tqdm(sample_list):
    print(f"  Sample: {sample_id}")
    anndata_path = anndata_dir / f"{sample_id}_with_interactions.pkl"

    with open(anndata_path, 'rb') as f:
        adata = pickle.load(f)

    if f"{overall_key}" in adata.uns:
        interaction_freqs = adata.uns[f"{overall_key}"]
        interaction_df_long = interaction_freqs.reset_index().melt(id_vars='index')
        interaction_df_long.columns = ['Cell_Type_1', 'Cell_Type_2', 'Interaction_Frequency']
        interaction_df_long['Sample_ID'] = sample_id
        # interaction_df_long['Niche_ID'] = niche_id
        #interaction_df_long.set_index(['Cell_Type_1', 'Cell_Type_2'], inplace=True) 
        overall_interactions.append(interaction_df_long)
    else:
        print(f"Interaction frequencies not found in sample {sample_id}")

#%%
overall_interaction_freqs = pd.concat(overall_interactions, ignore_index=True)
## replace NA values with 0, meaning no interactions observed for that pair in that sample
overall_interaction_freqs['Interaction_Frequency'] = overall_interaction_freqs['Interaction_Frequency'].fillna(0)
overall_medians = overall_interaction_freqs.groupby(['Cell_Type_1', 'Cell_Type_2'])['Interaction_Frequency'].median()
overall_interaction_freqs.set_index(['Cell_Type_1', 'Cell_Type_2'], inplace=True)


#%%
save_dir = output_dir / "redo" / "per_niche_lfc_above_median"
save_dir.mkdir(parents=True, exist_ok=True)
dict_results = {}
for niche_id in tqdm(niche_ids):
    print(f"Processing niche: {niche_id}")
    key_name = f'interaction_freqs_niche_{niche_id}'


    interactions = []
    above_median = pd.DataFrame(index=overall_medians.index, columns=['Above_Median'])
    above_median['Above_Median'] = 0
    below_median = pd.DataFrame(index=overall_medians.index, columns=['Below_Median'])
    below_median['Below_Median'] = 0
    num_samples_with_niche = 0
    for sample_id in sample_list:
        print(f"  Sample: {sample_id}")
        anndata_path = anndata_dir / f"{sample_id}_with_interactions.pkl"

        with open(anndata_path, 'rb') as f:
            adata = pickle.load(f)

        if f"{key_name}" in adata.uns:
            interaction_freqs = adata.uns[f"{key_name}"]
            interaction_df_long = interaction_freqs.reset_index().melt(id_vars='index')
            interaction_df_long.columns = ['Cell_Type_1', 'Cell_Type_2', 'Interaction_Frequency']
            # interaction_df_long['Sample_ID'] = sample_id
            # interaction_df_long['Niche_ID'] = niche_id
            interaction_df_long.set_index(['Cell_Type_1', 'Cell_Type_2'], inplace=True) 
            assert overall_medians.shape[0] == interaction_df_long.shape[0], "Mismatch in interaction frequencies shape."
            medians_aligned, interaction_df_long = overall_medians.align(interaction_df_long, join='inner', axis=0)
            interaction_df_long['lfc'] = np.log2((interaction_df_long['Interaction_Frequency'] + 1e-9) / (medians_aligned + 1e-9))
            # compare to overall median
            comparison = interaction_df_long['Interaction_Frequency'] > medians_aligned
            comparison_negative = interaction_df_long['Interaction_Frequency'] < medians_aligned
            above_median.loc[comparison.index, 'Above_Median'] += comparison.astype(int)
            below_median.loc[comparison_negative.index, 'Below_Median'] += comparison_negative.astype(int)
            interactions.append(interaction_df_long)
            num_samples_with_niche += 1
        else:
            print(f"Interaction frequencies for niche {niche_id} not found in sample {sample_id}")

    above_median['Above_Median_Fraction'] = above_median['Above_Median'] / num_samples_with_niche
    below_median['Below_Median_Fraction'] = below_median['Below_Median'] / num_samples_with_niche


    pseudocount = 1e-6
    mean_interactions = sum(interactions)/len(interactions)
    
    df_plot = mean_interactions.merge(above_median, left_index=True, right_index=True).merge(below_median, left_index=True, right_index=True)
    df_plot = df_plot.reset_index()
    ## if mean is greater than overall median, lfc is positive, else negative Interaction_Frequency_mean
    df_plot['size'] = np.where(df_plot['lfc'] > 0, df_plot['Above_Median_Fraction'], df_plot['Below_Median_Fraction'])
    dict_results[niche_id] = df_plot

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create scatter plot
    scatter = ax.scatter(
        x=df_plot['Cell_Type_2'],
        y=df_plot['Cell_Type_1'],
        s=df_plot['Above_Median_Fraction'] * 120,  # scale size for visibility
        c=df_plot['lfc'],
        cmap='coolwarm',
        alpha=0.8,
        edgecolor='k',
        linewidth=0.5,
        vmin=-3, vmax=3
    )

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Log2 Fold Change (LFC)', fontsize=12)

    # Titles and labels
    ax.set_title(f'Niche {niche_id}: Interactions (color = LFC, size = % support)', fontsize=14, pad=15)
    ax.set_xlabel('Cell Type (Column)', fontsize=12)
    ax.set_ylabel('Cell Type (Row)', fontsize=12)

    # Format ticks
    ax.set_xticks(range(len(df_plot['Cell_Type_2'].unique())))
    ax.set_yticks(range(len(df_plot['Cell_Type_1'].unique())))
    ax.set_xticklabels(df_plot['Cell_Type_2'].unique(), rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(df_plot['Cell_Type_1'].unique(), fontsize=8)

    # Grid and layout
    ax.grid(alpha=0.2, linestyle='--', linewidth=0.5)
    ax.set_facecolor('#f9f9f9')
    ax.set_aspect('equal', adjustable='box')  # keeps circles round and layout proportional
    plt.tight_layout()
    fig_path = save_dir / f"niche_{niche_id}_interaction_scatter.png"
    #plt.savefig(fig_path, dpi=300)
    plt.show()


# %%
results_dir = save_dir / "dataframes"
results_dir.mkdir(parents=True, exist_ok=True)
for niche_id, df_result in dict_results.items():
    df_path = results_dir / f"{niche_id}.parquet"
    df_result.to_parquet(df_path, engine='fastparquet')
# %%
cluster_path_old = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/annotation/clusters_annotated.parquet")
cluster_path_new = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/annotation/clusters_annotated_v2.parquet")

df_clusters_old = pd.read_parquet(cluster_path_old, engine='fastparquet')
df_clusters_new = pd.read_parquet(cluster_path_new, engine='fastparquet')

##mapping from old to new
#
# %%
orig_col = 'kmeans_24_seed686_filtered'
# unique pairs of niche col and original col
old_pairs = df_clusters_old[[orig_col, 'niche']].drop_duplicates()
new_pairs = df_clusters_new[[orig_col, 'niche']].drop_duplicates()
mapping = {}
for _, row in old_pairs.iterrows():
    orig_value = row[orig_col]
    niche_value = row['niche']
    new_niche_value = new_pairs[new_pairs[orig_col] == orig_value]['niche'].values
    if len(new_niche_value) > 0:
        mapping[niche_value] = new_niche_value[0]
    else:
        print(f"No match found for original value {orig_value} in new pairs.")
#%%
new_results_dir = save_dir / "dataframes_v2"
new_results_dir.mkdir(parents=True, exist_ok=True)
for niche_id, df_result in dict_results.items():
    new_niche_id = mapping.get(niche_id, None)
    if new_niche_id is not None:
        df_path = new_results_dir / f"{new_niche_id}.parquet"
        df_result.to_parquet(df_path, engine='fastparquet')
    else:
        print(f"No mapping found for niche ID {niche_id}, skipping save.")
# %%
