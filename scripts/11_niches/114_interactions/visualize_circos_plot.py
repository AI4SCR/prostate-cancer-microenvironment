#%%
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sys.path.append("/users/mensmeng/workspace/nhoods/PCa/05_nhoods/PCA_NHOODs_clean/05_niche_identification/visualization")
from utils_colors import color_dict_label, color_dict_niche
from circos_plots import plot_circos_plot


#%%
data_dir = Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/visualization/interactions/redo/per_niche_lfc_above_median/dataframes_v2")
#
niches = list(color_dict_niche.keys())
niche_colors = [color_dict_niche[niche] for niche in niches]
labels = list(color_dict_label.keys())
label_colors = [color_dict_label[label] for label in labels]


#%%
for niche in niches:
    df_path = data_dir / f"{niche}.parquet"
    df_example = pd.read_parquet(df_path, engine='fastparquet')
    # filter nas
    df_example = df_example.dropna()
    ##
    df_circos = df_example.set_index(['Cell_Type_1', 'Cell_Type_2'])

    ## rename columns  "Interaction_Frequency" to aggregated_interaction
    ##Above_Median_Fraction above_median_fraction
    df_circos.rename(columns={'Interaction_Frequency': 'aggregated_interaction', 'Above_Median_Fraction': 'above_median_fraction'}, inplace=True)
    df_circos = df_circos[['aggregated_interaction', 'above_median_fraction']]
    fig = plot_circos_plot(df=df_circos, color = 'above_median_fraction', color_map= color_dict_label, niche_name = niche, aggregator = None, save_path= None)
    plot_dir = data_dir.parent / "circos_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig_path = plot_dir / f"{niche}_circos_plot.pdf"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
# %%
