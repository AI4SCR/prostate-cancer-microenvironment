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

# Non-biological channels present in intensity tables: DNA intercalator +
# cell-segmentation-kit channels, plus FAP (excluded by the paper's Methods
# due to non-specific staining after in-house conjugation). Everything else
# in the panel is one of the 34 analysis markers.
NON_MARKER_CHANNELS = ["dna1", "dna2", "icsk1", "icsk2", "icsk3", "fap"]
# Non-marker columns carried alongside intensities by prepare_data()/export_for_r.py.
INDEX_COLUMNS = ["sample_id", "object_id", "slide_code", "donor_block_id", "pat_id"]
# Cell-type-label columns present in metadata.parquet (see ai4bmr_datasets.PCa.label_transfer()).
LABEL_COLUMNS = ["label", "main_group", "label_id", "main_group_id", "meta_label", "meta_label_id"]


def load_exported_cells(export_dir: Path, exclude_undefined: bool = False) -> pd.DataFrame:
    """Load and merge the per-cell tables `export_for_r.py` writes to `EXPORT_DIR`.

    One row per cell, combining `metadata.parquet` (labels) and
    `intensity_normalized.parquet` (markers), used by every figure script
    that plots or clusters individual cells.

    `exclude_undefined=True` drops the ~3% of cells the paper's Results
    describes as "remained unclassified and was excluded from the analysis"
    -- appropriate for cell-type-level figures (Figure 2a's heatmap), but
    NOT for whole-dataset figures like Figure 2b's UMAP, whose own caption
    states it covers all 2,191,967 cells. Decide per figure, don't default
    to filtering.
    """
    metadata = pd.read_parquet(export_dir / "metadata.parquet").reset_index()
    intensity = pd.read_parquet(export_dir / "intensity_normalized.parquet").reset_index()
    cells = intensity.merge(metadata, on=["sample_id", "object_id"], validate="one_to_one")
    if exclude_undefined:
        cells = cells[cells["label"] != "undefined"]
    return cells


def marker_columns(cells: pd.DataFrame, expected: int | None = 34) -> list[str]:
    """Return the marker-intensity columns of a `load_exported_cells()` table.

    `expected` asserts the resulting count (34 for the full panel); pass
    `None` to skip the check, e.g. when a figure restricts to a marker subset.
    """
    cols = [c for c in cells.columns if c not in NON_MARKER_CHANNELS + INDEX_COLUMNS + LABEL_COLUMNS]
    if expected is not None:
        assert len(cols) == expected, f"expected {expected} markers, got {len(cols)}: {cols}"
    return cols


def resolve_base_dir() -> Path:
    """Load `.env` and return `BASE_DIR` as a `Path`.

    Single source of truth for where scripts default to when `--base_dir` is
    not passed on the CLI, so every script points at the same dataset without
    a hardcoded machine-specific path. See REPRODUCIBILITY.md.
    """
    import os
    from dotenv import load_dotenv

    load_dotenv()
    base_dir = os.environ.get("BASE_DIR")
    assert base_dir, "BASE_DIR is not set; copy .env.example to .env and fill it in"
    return Path(base_dir).expanduser()


def resolve_export_dir() -> Path:
    """Load `.env` and return `EXPORT_DIR` as a `Path`, guaranteed outside `BASE_DIR`.

    `BASE_DIR` (the PCa dataset folder) is read-only for every script in this
    repo — nothing may ever write there. `EXPORT_DIR` is where generated
    tables/figures go instead. See REPRODUCIBILITY.md.
    """
    import os
    from dotenv import load_dotenv

    load_dotenv()
    export_dir = os.environ.get("EXPORT_DIR")
    assert export_dir, "EXPORT_DIR is not set; copy .env.example to .env and fill it in"
    return assert_outside_base_dir(Path(export_dir).expanduser())


def resolve_legacy_dir() -> Path:
    """Load `.env` and return `LEGACY_DATA_DIR` as a `Path`.

    Read-only, consolidated copy of pre-migration data with no reproducing
    script in this repo (niche-neighborhood raw data, the PCA_NHOODs_clean
    code some figure scripts import, colormaps.yaml). Never use this for
    metadata/clinical/intensity(_normalized).parquet -- those come from
    `resolve_export_dir()` instead. See REPRODUCIBILITY.md and
    figure_script_mapping.md.
    """
    import os
    from dotenv import load_dotenv

    load_dotenv()
    legacy_dir = os.environ.get("LEGACY_DATA_DIR")
    assert legacy_dir, "LEGACY_DATA_DIR is not set; copy .env.example to .env and fill it in"
    return Path(legacy_dir).expanduser()


def assert_outside_base_dir(path: Path) -> Path:
    """Fail loudly if `path` is inside `BASE_DIR` (or a symlink alias of it).

    `BASE_DIR` is the read-only PCa dataset folder on shared storage — no
    script in this repo may write there, directly or via a symlinked path
    (e.g. a home-directory shortcut into the same underlying tree).
    """
    import os

    base_dir = resolve_base_dir()
    real_path = os.path.realpath(path)
    real_base_dir = os.path.realpath(base_dir)
    assert not (real_path == real_base_dir or real_path.startswith(real_base_dir + os.sep)), (
        f"refusing to write inside BASE_DIR (the PCa dataset folder): {path} "
        f"resolves under {real_base_dir}. Use a path outside BASE_DIR."
    )
    return path


def get_colormap_dict(name: str, as_rgb: bool = False):
    import os
    import yaml
    from pathlib import Path

    colormaps_path = Path(__file__).resolve().parents[2] / "colormaps.yaml"
    if not colormaps_path.exists():
        base_dir = os.environ.get("BASE_DIR")
        assert base_dir, "colormaps.yaml not found next to the repo; set BASE_DIR or add colormaps.yaml at the repo root"
        colormaps_path = Path(base_dir) / "colormaps.yaml"

    with open(colormaps_path) as f:
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


def prepare_data(base_dir: Path, scale="minmax", mask_version: str = "annotated"):
    """Load and normalize cell intensities from the PCa dataset.

    `mask_version="filtered"` gives every segmented cell (no cell-type labels
    yet) and is the correct input for the 01-clustering scripts that assign
    labels for the first time. `mask_version="annotated"` (default) requires
    `01_raw/annotations/labels.parquet` to already exist and is what every
    downstream figure script should use once labels have been produced.
    See REPRODUCIBILITY.md for the full bootstrap order.
    """
    from ai4bmr_datasets import PCa

    dataset = PCa(base_dir=base_dir,
             image_version='filtered',
             mask_version=mask_version,
             load_intensity=True,
             load_metadata=False,  # only .intensity and .clinical are used below
             align=False)

    dataset.setup()

    data = dataset.intensity
    assert data.isna().any().any() == False
    metadata = dataset.clinical

    metadata_cols = ["slide_code", "donor_block_id", "pat_id"]
    assert metadata[metadata_cols].isna().sum().sum() == 0

    data = data.join(metadata[metadata_cols]).set_index(metadata_cols, append=True)
    data = data.sort_index(level=["sample_id"])

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

def normalize_img(img, censoring=0.999, cofactor=1, exclude_zeros=True):
    """Normalizes an image using an arcsin/work/FAC/FBM/DBC/mrapsoma/prometex/data/omics-embed/datasetsh transformation and applies intensity censoring.

    Args:
        img (np.ndarray): The input image array.
        censoring (float, optional): The quantile to censor the image intensities. Defaults to 0.999.
        cofactor (int, optional): The cofactor for the arcsinh transformation. Defaults to 1.
        exclude_zeros (bool, optional): Whether to exclude zero values when computing the censoring threshold. Defaults to True.

    Returns:
        np.ndarray: The normalized and censored image array.
    """
    img = np.arcsinh(img / cofactor)

    if exclude_zeros:
        masked_img = np.where(img == 0, np.nan, img)
        thres = np.nanquantile(masked_img, censoring, axis=(1, 2), keepdims=True)
        # All-zero channels become all-NaN after masking; use 0 threshold fallback.
        thres = np.where(np.isnan(thres), 0.0, thres)
    else:
        thres = np.quantile(img, q=censoring, axis=(1, 2), keepdims=True)

    img = np.minimum(img, thres)

    return img