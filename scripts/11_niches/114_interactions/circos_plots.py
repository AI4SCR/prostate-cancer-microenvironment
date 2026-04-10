#%%
from pycirclize import Circos
from matplotlib import cm, colors
from typing import Dict
import pandas as pd
import numpy as np

#%%
def check_triangular_zero_xor(df: pd.DataFrame):
    """
    Checks if the strictly upper triangle OR the strictly lower triangle of a 
    square DataFrame is composed entirely of zeros, but NOT both.

    Args:
        df: The input Pandas DataFrame (assumed to be square).

    Returns:
        Boolean result for check.
    """
    if df.shape[0] != df.shape[1]:
        return False, "DataFrame is not square. Cannot perform triangular check."

    # Convert the DataFrame to a NumPy array for efficient checks
    data = df.values

    # Upper Triangle 
    is_upper_zero = np.all(np.triu(data, k=1) == 0) #k=1 to exclude the diagonal itself.
    # Lower Triangle
    is_lower_zero = np.all(np.tril(data, k=-1) == 0) #k=-1 to exclude the diagonal itself.

    # 3. Apply the XOR (Exclusive OR) logic: A OR B, but NOT (A AND B)
    result_xor = is_upper_zero ^ is_lower_zero

    return result_xor


def get_color_map(df: pd.DataFrame) -> Dict[str, str]:
    '''
    Generate a color map for cell types in the dataframe.
    Args:
        df : pandas.DataFrame -> dataframe containing input data.
        It has to contain:
                - MultiIndex with 'cell_type_1' and 'cell_type_2'(str).
    returns:
        Dict[str, str]: Dictionary mapping cell types to colors.
    '''
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    if df.index.names != ['cell_type_1', 'cell_type_2']:
        df.index.set_names(['cell_type_1', 'cell_type_2'], inplace=True)

    df_new = df.reset_index()
    labels = list(set(df_new['cell_type_1'].unique()).union(set(df_new['cell_type_2'].unique())))
    cmap = cm.get_cmap('tab10')
    color_indices = {label: i for i, label in enumerate(labels)}
    color_map = {ct: cmap(color_indices[ct]) for ct in labels}
    return color_map


def data_for_circos_plot(df: pd.DataFrame, color: str):
    '''
    Create color and width dictionaries for circos plot links. 
        they will contain:
            labels as keys (tuples) and color/width values as values.    
    Creates sectors dataframe with width values.
        It will contain:
            labels as index and columns, and width values as values.
            The dataframe has 0 values on the top (or bottom) of the matrix -> not count interactions twice.

    Args:
        df : pandas.DataFrame -> dataframe containing input data.
            It has to contain:
                - MultiIndex with 'cell_type_1' and 'cell_type_2'(str).
                - column 'above_median_fraction' (float): Fraction of samples with interaction above median.
                - column 'aggregated_interaction' (float): Aggregated interaction value.
        color: String indicating which color values to use 'aggregated_interaction' or 'above_median_fraction'
    Returns:
        Dict[str, Dict[Tuple[str, str], float]]: Dictionary with 'color_dict' and 'width_dict'.
    '''
    assert color in ['aggregated_interaction', 'above_median_fraction'], "color must be either 'aggregated_interaction' or 'above_median_fraction'"
    assert 'above_median_fraction' and 'aggregated_interaction' in df.columns, "DataFrame must contain 'above_median_fraction' and 'aggregated_interaction' columns"
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    if df.index.names != ['cell_type_1', 'cell_type_2']:
        df.index.set_names(['cell_type_1', 'cell_type_2'], inplace=True)

    #WIDTH
    width = [col for col in df.columns.tolist() if col != color][0]
    width_df = df[width]
    width_df = width_df.reset_index()
    width_df.columns = ['from', 'to', 'Value']
    width_df['Value'].dropna()
    width_dict = {(row['from'], row['to']): row['Value'] for _, row in width_df.iterrows()}

    # COLOR
    color_df = df[color]
    color_df = color_df.reset_index()
    color_df.columns = ['from', 'to', 'Value']
    color_df['Value'].dropna()
    color_dict = {(row['from'], row['to']): row['Value'] for _, row in color_df.iterrows()}

    #SECTORS 
    sectors_df = df[width].copy()
    sectors_df = sectors_df.reset_index()
    sectors_df.columns = ['cell_type_1', 'cell_type_2', 'width']
    sectors_df = sectors_df.pivot(index='cell_type_2', columns='cell_type_1', values='width')
    sectors_df = sectors_df.fillna(0)
    assert check_triangular_zero_xor(sectors_df), "Either the upper or lower triangle values of sectors_df must be 0, but not both."

    return {'color_dict': color_dict, 'width_dict': width_dict, 'sectors_df': sectors_df}


def get_circos(df: pd.DataFrame, color:str = 'aggregated_interaction', color_map: Dict[str,str] = None):
    '''
    Create circos object for aggregated interactions across samples.

    Args:
        df : pandas.DataFrame -> dataframe containing input data.
            It has to contain:
                - MultiIndex with 'cell_type_1' and 'cell_type_2'(str).
                - column 'above_median_fraction' (float): Fraction of samples with interaction above median.
                - column 'aggregated_interaction' (float): Aggregated interaction value.
        color: String indicating which color values to use 'aggregated_interaction' or 'above_median_fraction'
        color_map: Dictionary mapping cell types to colors. Default is None, a default colormap is generated.
    Returns:
        Circos object
    '''
    assert color in ['aggregated_interaction', 'above_median_fraction'], "color must be either 'aggregated_interaction' or 'above_median_fraction'"
    assert 'above_median_fraction' and 'aggregated_interaction' in df.columns, "DataFrame must contain 'above_median_fraction' and 'aggregated_interaction' columns"
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    dicts = data_for_circos_plot(df = df, color = color)

    color_dict = dicts['color_dict']
    width_dict = dicts['width_dict']
    sectors_df = dicts['sectors_df']
    
    if color_map is None:
        cell_color_map= get_color_map(df)
    else:
        cell_color_map = color_map
    
    color_var = color
    val_min = min(color_dict.values())
    val_max = max(color_dict.values())

    # define color and with of the links
    def link_handler(from_label, to_label):
        # Get external value
        val = color_dict.get((from_label, to_label)) or color_dict.get((to_label, from_label), 0)
        lw = width_dict.get((from_label, to_label)) or width_dict.get((to_label, from_label), 1)

        if color_var == 'aggregated_interaction':
            filt = val
        else:
            filt = lw
        
        if filt <= 0.001:
            # Filter = if the value used for the color (either aggregated interaction or above median fraction) is 0.001 or less, make it invisible
            color = 'none'
            lw = 0
        else:    
            # If the value used for the color is bigger than 0.001 => assign it to a color using a colormap
            cmap = cm.get_cmap("Reds")  # or "Reds", "coolwarm", etc.
            val_min = min(color_dict.values())
            val_max = max(color_dict.values())
            norm = colors.Normalize(vmin=val_min, vmax=val_max)
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            color = sm.to_rgba(val)  # val should be in 0-1

        # Return styling dictionary
        return dict(ec='none', lw=lw, fc=color, alpha=0.7)   
    
    circos = Circos.chord_diagram(
        sectors_df,
        space=2,
        cmap = cell_color_map,
        label_kws=dict(
            size=14,          # larger font
            color="black",    # dark text
            r=110,            # radial distance from circle
            orientation="vertical",  # or 'horizontal'
        ),
        link_kws_handler=link_handler
        )
    return circos
   

def plot_circos_plot(df: pd.DataFrame, color = 'aggregated_interaction', color_map: Dict[str,str] = None, niche_name: str = None, aggregator: str = None, save_path: str = None):
    '''
    Create circos plots for aggregated interactions across samples.

    Args:
        df : pandas.DataFrame -> dataframe containing input data.
            It has to contain:
                - MultiIndex with 'cell_type_1' and 'cell_type_2'(str).
                - column 'above_median_fraction' (float): Fraction of samples with interaction above median.
                - column 'aggregated_interaction' (float): Aggregated interaction value.
        color: String indicating which color values to use 'aggregated_interaction' or 'above_median_fraction'
        color_map: Dictionary mapping cell types to colors. Default is None, a default colormap is generated.
        niche_name: String indicating the name of the niche to put in the title if wanted. Default is None.
        aggregator: String indicating the aggregation method used to put in the title if wanted. Default is None.
        save_path: String indicating the path to save the figure. If None, figure is not saved. Default is None.
    Returns:
        None: Displays dot plot.
        Saves figure if save_path is provided.
    '''
    assert color in ['aggregated_interaction', 'above_median_fraction'], "color must be either 'aggregated_interaction' or 'above_median_fraction'"
    assert 'above_median_fraction' and 'aggregated_interaction' in df.columns, "DataFrame must contain 'above_median_fraction' and 'aggregated_interaction' columns"
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    circos= get_circos(df=df, color=color, color_map=color_map) #circo object

    fig = circos.plotfig(figsize=(20, 15))

    # Adjust spacing to make room for title and colorbar
    fig.subplots_adjust(top=3, right=3)  # leave margin on top and right
    # Title 
    if niche_name==None:
        niche_name=''
    if aggregator==None:
        aggregator=''
    else:
        aggregator=f'({aggregator}) '
    fig.suptitle(
        f"Circos Plot of Aggregated Interactions {aggregator}in Niche {niche_name}\nColor = {color}, Line Width = {[col for col in df.columns.tolist() if col != color][0]}",
        fontsize=16,
        fontweight="bold",
        y=0.99  # vertical position: 1.0 is top of figure
    )
    

    # Create colorbar
    val_min = min(df[color])
    val_max = max(df[color])
    norm = colors.Normalize(vmin=val_min, vmax=val_max)
    cmap = cm.get_cmap("Reds")
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar_ax = fig.add_axes([0.91, 0.25, 0.02, 0.5])  # [left, bottom, width, height]
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical', label=color)
    cbar.set_label(color, fontsize=14)
    cbar.ax.tick_params(labelsize=12)
    fig.show()    

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig