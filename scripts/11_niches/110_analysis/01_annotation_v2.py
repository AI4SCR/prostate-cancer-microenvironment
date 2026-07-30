#%%
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

#%%
output_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/annotation")
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Results will be saved to {output_dir}")


# #%%
# annotation_dict_niche = {
#     "0": 'luminal',
#     "1": 'luminal_CAF2s_CAF1(CD105+)',
#     "2": 'tumor(ERG+)',
#     "3": 'endo_luminal',
#     "4": 'CAFs_CAF2(AR-)_CAF1(CD105-)_fibrocytes_Tcells',
#     "5": 'canonical_BLepithelium',
#     "6": 'peritumoral_CAFs_CAF1(CD105+)_Tcells_endo',
#     "7": 'luminal',
#     "8": 'TLS_Bcells_Tcells',
#     "9": 'tumorERG+p53+_ProlifLuminal',
#     "10": 'CAF2s_CAF2(AR+)_CAF1(CD105-)_Tcells',
#     "11": 'endo_luminal',
#     "12": 'luminal_CAF1(CD105-)_undefined',
#     "13": 'CAF2s_CAF2(AR+CES1+)_CAF2(AR+)',
#     "14": 'endo_Tcells_CAF2s_CAF1(CD105-)',
#     "16": 'Macrophages_Tcells_CAF2s_CAF1(CD105-)',
#     "17": 'peritumoral_CAF1(CD105+)',
#     "18": 'luminal_lumERG+_WellDiffTumors',
#     "19": 'canonical_BLepithelium',
#     "20": 'vessels_pericytes_CAF1(CD105-)',
#     "21": 'CAF2s_CAF2(AR-)_CAF1(CD105-)',
#     "22": 'endo_Tcells_CAF2s_CAF1(CD105-)',
#     "23": 'CAF2s_CAF2(AR-)_CAF1(CD105-)',
#     "unassigned": 'unassigned'
# }
#%%
cluster_path = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/clustering/clusters.parquet")
df_clusters = pd.read_parquet(cluster_path, engine='fastparquet')
print(f"Clusters shape: {df_clusters.shape}")

#%% read dicts from excel table
annot_excel = Path("/users/mensmeng/workspace/nhoods/PCa/05_nhoods/PCA_NHOODs_clean/niche_annotations_revised.xlsx")
df_anno = pd.read_excel(annot_excel)

annotation_dict_niche = dict(zip(df_anno['cluster'].astype(str), df_anno['niche']))
annotation_dict_meta_niche = dict(zip(df_anno['niche'], df_anno['meta_niche']))
color_dict_niche = dict(zip(df_anno['niche'], df_anno['niche_color']))
color_dict_meta_niche = dict(zip(df_anno['meta_niche'], df_anno['meta_niche_color']))

# %%
cluster_name = df_clusters.columns[0]
df_clusters[cluster_name] = df_clusters[cluster_name].astype(str)
df_clusters['niche'] = df_clusters[cluster_name].map(annotation_dict_niche).fillna('unassigned')
# # %%
# annotation_dict_meta_niche = {
#     'luminal': 'luminal_other_epithelial_cells',
#     'luminal_CAF2s_CAF1(CD105+)': 'luminal_CAFs',
#     'tumor(ERG+)': 'tumor_cells_+/-_CAFs',
#     'endo_luminal': 'luminal_other_epithelial_cells',
#     'CAFs_CAF2(AR-)_CAF1(CD105-)_fibrocytes_Tcells': 'CAF_subtypes_+/-_immune',
#     'canonical_BLepithelium': 'luminal_other_epithelial_cells',
#     'peritumoral_CAFs_CAF1(CD105+)_Tcells_endo': 'tumor_cells_+/-_CAFs',
#     'TLS_Bcells_Tcells': 'immune_cells',
#     'tumorERG+p53+_ProlifLuminal': 'tumor_cells_+/-_CAFs',
#     'CAF2s_CAF2(AR+)_CAF1(CD105-)_Tcells': 'CAF_subtypes_+/-_immune',
#     'luminal_CAF1(CD105-)_undefined': 'luminal_CAFs',
#     'CAF2s_CAF2(AR+CES1+)_CAF2(AR+)': 'CAF_subtypes_+/-_immune',
#     'endo_Tcells_CAF2s_CAF1(CD105-)': 'CAF_subtypes_+/-_immune',
#     'Macrophages_Tcells_CAF2s_CAF1(CD105-)': 'immune_cells',
#     'peritumoral_CAF1(CD105+)': 'tumor_cells_+/-_CAFs',
#     'luminal_lumERG+_WellDiffTumors': 'tumor_cells_+/-_CAFs',
#     'vessels_pericytes_CAF1(CD105-)': 'luminal_CAFs',
#     'CAF2s_CAF2(AR-)_CAF1(CD105-)': 'CAF_subtypes_+/-_immune',
#     'CAF2s_CAF2(AR-)_CAF1(CD105-)': 'CAF_subtypes_+/-_immune',
#     'unassigned': 'unassigned'
# }

# # %%
# color_dict_meta_niche = {
#     'luminal_other_epithelial_cells': "#e12626",
#     'luminal_CAFs': "#fd79f4",
#     'tumor_cells_+/-_CAFs': "#9905d4",
#     'CAF_subtypes_+/-_immune': '#1f77b4',
#     'immune_cells': '#2ca02c',
#     'unassigned': '#7f7f7f'
# }


# # %%
# color_dict_niche = {
#     'luminal': '#e12626',
#     'luminal_CAF2s_CAF1(CD105+)': "#f6b4f2",
#     'tumor(ERG+)': "#360549",
#     'endo_luminal': "#ea7e7e",
#     'CAFs_CAF2(AR-)_CAF1(CD105-)_fibrocytes_Tcells': "#448cbf",
#     'canonical_BLepithelium': "#e19f26",
#     'peritumoral_CAFs_CAF1(CD105+)_Tcells_endo': '#9905d4',
#     'TLS_Bcells_Tcells': "#22985D",
#     'tumorERG+p53+_ProlifLuminal': "#91617b",
#     'CAF2s_CAF2(AR+)_CAF1(CD105-)_Tcells': "#7d8ec2",
#     'luminal_CAF1(CD105-)_undefined': '#fd79f4',
#     'CAF2s_CAF2(AR+CES1+)_CAF2(AR+)': "#6ed5e9",
#     'endo_Tcells_CAF2s_CAF1(CD105-)': "#530ce3",
#     'Macrophages_Tcells_CAF2s_CAF1(CD105-)': "#1ada1a",
#     'peritumoral_CAF1(CD105+)': "#763457",
#     'luminal_lumERG+_WellDiffTumors': "#d44005",
#     'vessels_pericytes_CAF1(CD105-)': "#AA7804",
#     'CAF2s_CAF2(AR-)_CAF1(CD105-)': "#77abd1",
#     'unassigned': '#7f7f7f'
# }

# %%
## create dataframe out of annotation dict_niche with key as cluster and value as niche
df_annotation_niche = pd.DataFrame.from_dict(annotation_dict_niche, orient='index', columns=['niche'])
df_annotation_niche.index.name = 'cluster'
df_annotation_niche.reset_index(inplace=True)

## map meta_niche
df_annotation_niche['meta_niche'] = df_annotation_niche['niche'].map(annotation_dict_meta_niche).fillna('unassigned')

## map colors
df_annotation_niche['niche_color'] = df_annotation_niche['niche'].map(color_dict_niche).fillna('#7f7f7f')
df_annotation_niche['meta_niche_color'] = df_annotation_niche['meta_niche'].map(color_dict_meta_niche).fillna('#7f7f7f')
# %%
df_annotation_niche.to_csv(output_dir / 'niche_annotations_v2.csv', index=False)

# # %%
# %%
cluster_name = df_clusters.columns[0]
df_clusters[cluster_name] = df_clusters[cluster_name].astype(str)
df_clusters['niche'] = df_clusters[cluster_name].map(annotation_dict_niche).fillna('unassigned')

df_clusters['meta_niche'] = df_clusters['niche'].map(annotation_dict_meta_niche).fillna('unassigned')

#%%
df_clusters['niche'] = df_clusters['niche'].astype('category')
df_clusters['meta_niche'] = df_clusters['meta_niche'].astype('category')
# %%
df_clusters.to_parquet(output_dir / 'clusters_annotated_v2.parquet', engine='fastparquet', index=True)
# %%
