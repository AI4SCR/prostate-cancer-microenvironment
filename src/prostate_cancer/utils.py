from pathlib import Path

import colorcet as cc
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_rgba
from pandas.api.types import is_numeric_dtype
from sklearn.preprocessing import MinMaxScaler, StandardScaler

def get_colormap_dict(name: str, as_rgb: bool = False):
    import yaml
    with open('/work/FAC/FBM/DBC/mrapsoma/prometex/projects/PCa/colormaps.yaml') as f:
        colormaps = yaml.load(f, Loader=yaml.SafeLoader)
    colormap = colormaps[name]
    return colormap

def normalize_data(data, scale="minmax"):
    # NORMALIZE

    censoring = 0.999
    cofactor = 1

    x = np.arcsinh(data / cofactor)
    thres = np.quantile(x, censoring, axis=0)

    for idx, t in enumerate(thres):
        x.values[:, idx] = np.where(x.values[:, idx] > t, t, x.values[:, idx])

    if scale == "minmax":
        data = pd.DataFrame(
            MinMaxScaler().fit_transform(x), columns=x.columns, index=x.index
        )
    elif scale == "standard":
        data = pd.DataFrame(
            StandardScaler().fit_transform(x), columns=x.columns, index=x.index
        )
    elif scale is None:
        data = x
    else:
        raise NotImplementedError()

    return data

def normalize(data: pd.DataFrame, scale: str = 'minmax', exclude_zeros: bool = False):
    import numpy as np

    index = data.index
    columns = data.columns
    x = data.values

    censoring = 0.999
    cofactor = 1
    x = np.arcsinh(x / cofactor)

    if exclude_zeros:
        masked_x = np.where(x == 0, np.nan, x)
        thres = np.nanquantile(masked_x, censoring, axis=0)
    else:
        thres = np.nanquantile(x, censoring, axis=0)

    x = np.minimum(x, thres)
    assert (x.max(axis=0) <= thres).all()

    if scale == "minmax":
        x = MinMaxScaler().fit_transform(x)
    elif scale == "standard":
        x = StandardScaler().fit_transform(x)
    else:
        raise NotImplementedError()

    return pd.DataFrame(x, index=index, columns=columns)


def prepare_data(base_dir: Path, scale="minmax"):
    from ai4bmr_datasets import PCa

    dataset = PCa(base_dir=base_dir,
             image_version='filtered',
             mask_version='annotated',
             load_intensity=True,
             load_metadata=True,
             align=False)

    dataset.setup()

    data = dataset.intensity
    assert data.isna().any().any() == False
    metadata = dataset.clinical

    metadata_cols = ["slide_code", "donor_block_id", "pat_id"]
    assert metadata[metadata_cols].isna().sum().sum() == 0

    data = data.join(metadata[metadata_cols]).set_index(metadata_cols, append=True)
    data = data.sort_index(level=["sample_name"])

    # NORMALIZE
    censoring = 0.999
    cofactor = 1
    x = np.arcsinh(data / cofactor)
    thres = np.quantile(x, censoring, axis=0)
    for idx, t in enumerate(thres):
        x.values[:, idx] = np.where(x.values[:, idx] > t, t, x.values[:, idx])

    if scale == "minmax":
        data = pd.DataFrame(
            MinMaxScaler().fit_transform(x), columns=x.columns, index=x.index
        )
    elif scale == "standard":
        data = pd.DataFrame(
            StandardScaler().fit_transform(x), columns=x.columns, index=x.index
        )
    elif scale is None:
        data = x
    else:
        raise NotImplementedError()

    return data


def create_color_maps_from_index(data: pd.DataFrame):
    color_maps = {}
    for label_name in set(data.index.names) - {"object_id"}:
        n = len(cc.glasbey_category10)
        labels = np.sort(data.index.get_level_values(label_name).unique())
        color_map = {
            label: cc.glasbey_category10[i % n] for i, label in enumerate(labels)
        }
        color_maps[label_name] = color_map
    return color_maps


def create_color_maps_from_frame(data: pd.DataFrame):
    color_maps = {}
    for col_name in set(data.columns):
        n = len(cc.glasbey_category10)
        labels = np.sort(data[col_name].unique())
        color_map = {
            label: cc.glasbey_category10[i % n] for i, label in enumerate(labels)
        }
        color_maps[col_name] = color_map
    return color_maps


def create_color_maps(data: pd.DataFrame):
    color_maps = {}
    for col in data:
        if data[col].dtype.name == "category":
            n = len(cc.glasbey_category10)
            labels = np.sort(data[col].unique())
            color_map = {
                label: to_rgba(cc.glasbey_category10[i % n])
                for i, label in enumerate(labels)
            }
            color_maps[col] = color_map
        elif is_numeric_dtype(data[col]):
            cmap = LinearSegmentedColormap.from_list(col, cc.linear_kry_0_97_c73)
            color_maps[col] = cmap
        else:
            print(f"WARNING: {col} is not a category or numeric dtype")
    return color_maps


def plot_umap_index(data: pd.DataFrame, embedding, color_maps: dict, save_dir: Path):
    label_names = set(data.index.names) - {"object_id"}
    for label_name in label_names:
        labels = data.index.get_level_values(label_name)

        color_map = color_maps[label_name]

        # create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(
            x=embedding[:, 0],
            y=embedding[:, 1],
            c=[color_map[l] for l in labels],
            s=1,
            alpha=0.3,
        )
        ax.set_facecolor("black")
        ax.set_title(label_name)

        # ax = umap.plot.points(mapper, labels=labels, background='black', show_legend=False)

        # create legend
        legend_patches = [
            mpatches.Patch(color=color, label=label)
            for label, color in color_map.items()
        ]
        ax.legend(handles=legend_patches, loc="center left", bbox_to_anchor=(1, 0.5))

        ax.figure.tight_layout()
        # ax.figure.show()
        ax.figure.savefig(save_dir / f"umap_{label_name}.png", dpi=300)
        plt.close(ax.figure)


def plot_umap_columns(data, embedding, save_dir: Path):
    value_names = set(data.columns)
    for value_name in value_names:
        values = data[value_name]
        cmap = matplotlib.colormaps.get_cmap("hot")

        # create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(x=embedding[:, 0], y=embedding[:, 1], c=cmap(values), s=1, alpha=0.3)
        ax.set_facecolor("black")

        ax.set_title(value_name)
        ax.figure.tight_layout()
        # ax.figure.show()
        ax.figure.savefig(save_dir / f"umap_{value_name}.png", dpi=300)
        plt.close(ax.figure)


def create_legends(row_annotations, color_maps):
    legends = []
    for i in row_annotations:
        type_ = "discrete" if row_annotations[i].dtype == "category" else "continuous"
        if type_ == "discrete":
            legend = {"type": type_, "label_to_color": color_maps[i], "title": i}
        else:
            legend = {
                "type": type_,
                "height": 75,
                "width": 10,
                "vmin": row_annotations[i].min(),
                "vmax": row_annotations[i].max(),
                "colormap": color_maps[i],
                "orientation": "vertical",
                "title": i,
            }
        legends.append(legend)

    return legends


def map_row_annotations_to_colors(row_data, color_maps):
    row_colors = row_data.copy()
    cat_cols = row_colors.select_dtypes("category").columns
    for label_name in cat_cols:
        cmap = color_maps[label_name]
        row_colors[label_name] = [cmap[v] for v in row_colors[label_name]]

    num_cols = row_data.select_dtypes(["float", "int"]).columns
    for label_name in num_cols:
        cmap = color_maps[label_name]
        row_colors[label_name] = [cmap(v) for v in row_colors[label_name]]
    return row_colors


def normalize_row_annotations(row_annotations):
    row_annotations_norm = row_annotations.copy()

    num_cols = row_annotations.select_dtypes(["float", "int"]).columns
    for col in num_cols:
        from matplotlib.colors import Normalize

        norm = Normalize(
            vmin=row_annotations[col].min(), vmax=row_annotations[col].max()
        )
        row_annotations_norm[col] = norm(row_annotations[col])
    return row_annotations_norm

