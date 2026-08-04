# %%
"""Reproduce Supplementary Figure 3a: TMA (core)-level cell-type composition
stacked barplot -- the `group_var == 'tma_id'` sibling of Figure 4a's
`pat_id` branch.

1:1 port of the old repo's
`000_paper/sync_paper/05-heterogeneity/stacked-frequencies-label.py`'s
`group_var == 'tma_id'` branch only -- not the `pat_id` branch (already
`figure4_patient_clustering.py`) and not the tail metagroup-averaging block
(only reachable for `pat_id`, already `figure4_metagroup_barplot.py`).
Duplicates the same helper functions
(`get_order`/`add_col_annotations`/`add_dendrogram_top`/
`get_label_frequency_table`) already in `figure4_patient_clustering.py`
rather than importing them, per this repo's own anti-pattern rule (no
cross-script imports of reusable logic) -- both copies are intentionally
identical.

Computes fresh (does not read the precomputed `LEGACY_DATA_DIR` ground
truth), matching Figure 4a's established pattern; nothing downstream reads
this script's own dendrogram-leaf-color output.

Reads from `$EXPORT_DIR` and `resources/colormaps.yaml`. Writes to
`$OUTPUT_FIGURES_DIR/figureS3/`.
"""
from pathlib import Path

import colorcet as cc
import matplotlib
import numpy as np
import pandas as pd
import yaml
from jsonargparse import CLI
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist

from prostate_cancer.utils import resolve_export_dir, resolve_output_figures_dir

matplotlib.use("Agg")

DISTANCE_THRESHOLD = 0.4  # height cut, per stacked-frequencies-label.py


def compute_label_frequency(data: pd.DataFrame, level: str, pseudocount: int = 1, group_vars: list[str] = ["sample_id"]) -> pd.Series:
    """1:1 port of datamodules/utils.py:compute_label_frequency()."""
    data[level] = data[level].astype("category")
    if pseudocount > 0:
        pdat = data.groupby(group_vars, observed=False)[level].value_counts()
        pdat += 1
        pdat /= pdat.groupby(group_vars, observed=False).sum()
        pdat.name = "proportion"
    else:
        pdat = data.groupby(group_vars, observed=False)[level].value_counts(normalize=True)
    return pdat


def get_label_frequency_table(data: pd.DataFrame, level: str, group_vars: list[str] = ["sample_id"]) -> pd.DataFrame:
    """1:1 port of datamodules/utils.py:get_label_frequency_table()."""
    props = compute_label_frequency(data=data, level=level, pseudocount=1, group_vars=group_vars)
    props = props.reset_index().pivot(index=group_vars, columns=level, values="proportion")
    props.columns = props.columns.astype(str)
    return props.astype(float)


def get_order(proportions: pd.DataFrame, method: str = "average"):
    dcond = pdist(proportions.values, metric="jensenshannon")
    Z = linkage(dcond, method=method)
    order_idx = leaves_list(Z)
    order = proportions.index[order_idx].tolist()
    return order, Z


def add_col_annotations(*, ax, mdat: pd.DataFrame, col_names: list[str], y0: float = 1.01, height: float = 0.03, colormaps=None):
    x0, y0, width, height = 0, y0, 1, height

    colors = []
    for i, col_name in enumerate(col_names):
        if col_name in colormaps:
            colormap_dict = colormaps[col_name]
            colormap_dict = {k: to_rgb(v) for k, v in colormap_dict.items()}
        else:
            n = len(cc.glasbey_bw)
            colormap_dict = {k: to_rgb(cc.glasbey_bw[i % n]) for i, k in enumerate(mdat[col_name].unique())}

        labels = mdat.reset_index()[col_name]
        colors.append(np.array([colormap_dict[str(label)] for label in labels]))

    colors = np.stack(colors)
    ax2 = ax.inset_axes([x0, y0, width, height])
    ax2.imshow(colors, aspect="auto")
    ax2.set_xticks([])
    ax2.set_yticks(range(len(col_names)))
    ax2.set_yticklabels(col_names)
    ax2.set_clip_on(False)
    return ax2


def add_dendrogram_top(ax, Z, thres: float = 0, y0: float = 1.01, height: float = 0.15):
    """Add a dendrogram aligned with the bar categories at the top of ax."""
    dend_ax = ax.inset_axes([0.0, y0, 1.0, height])
    dendro = dendrogram(
        Z,
        orientation="top",
        no_labels=False,
        color_threshold=thres,
        above_threshold_color="black",
        ax=dend_ax,
    )
    dend_ax.set_xticks([])
    dend_ax.set_xticklabels([])
    return dend_ax, dendro


def main(export_dir: Path | None = None):
    export_dir = export_dir or resolve_export_dir()
    save_dir = resolve_output_figures_dir() / "figureS3"
    save_dir.mkdir(parents=True, exist_ok=True)
    resources_dir = Path(__file__).resolve().parents[2] / "resources"

    logger.info("loading exported tables")
    df_labels = pd.read_parquet(export_dir / "metadata.parquet").reset_index()
    df_metadata = pd.read_parquet(export_dir / "clinical.parquet")

    # %% filter to tumor samples only
    df_sample_id = df_metadata.reset_index()
    df_sample_id = df_sample_id[~df_sample_id["is_tumor"].isna()]
    df_sample_id = df_sample_id[df_sample_id["is_tumor"] == "yes"]
    df_labels = df_labels.merge(df_sample_id, on="sample_id", how="inner")
    df_labels = df_labels.set_index(["sample_id", "object_id"])

    var_name = "label"
    group_var = "tma_id"
    tma_cols = ["stromogenic_smc_loss_reactive_stroma_present", "inflammation", "glandular_atrophy_pin", "gleason_grp", "pat_id"]

    valid_tma_ids = df_labels["tma_id"].unique().tolist()
    df_metadata = df_metadata[df_metadata["tma_id"].isin(valid_tma_ids)]
    df_metadata = df_metadata.set_index("tma_id")

    metadata = df_metadata.copy()
    metadata = metadata[tma_cols]
    metadata = metadata.reset_index()
    metadata = metadata.drop_duplicates(subset=[group_var])
    metadata = metadata.set_index(group_var)

    df_freqs = get_label_frequency_table(data=df_labels, level=var_name, group_vars=[group_var])
    df_freqs, metadata = df_freqs.align(metadata, join="inner", axis=0)
    logger.info(f"frequencies shape after alignment: {df_freqs.shape}")
    logger.info(f"metadata shape after alignment: {metadata.shape}")

    # %% cluster samples by JS distance of frequencies and add cluster annotation to metadata
    thresh = DISTANCE_THRESHOLD
    order, Z = get_order(df_freqs, method="average")
    fixed_clusters = fcluster(Z, t=thresh, criterion="distance")

    df_freqs, metadata = df_freqs.align(metadata, join="inner", axis=0)
    metadata["cluster"] = fixed_clusters  # computed, unused downstream -- kept verbatim per legacy source

    cols_to_plot = ["gleason_grp", "pat_id"]
    mdat = metadata.loc[order]
    props = df_freqs.loc[order]
    assert len(mdat) == len(props)

    # %% load colormaps
    with open(resources_dir / "colormaps.yaml") as f:
        colormaps = yaml.safe_load(f)
    colordict_level = colormaps[var_name]

    # %% stacked barplot of frequencies ordered by JS distance, with annotation bars above
    ax = props.plot(
        kind="bar",
        stacked=True,
        figsize=(25, 10),
        color=[colordict_level[c] for c in props.columns],
        width=0.8,
    )
    ax.set_ylim(0, 1.01)
    ax.set_ylabel("Proportion")
    ax.set_xlabel("TMA ID")
    ax.set_xticklabels("")
    ax.get_legend().remove()

    y0, height = 1.01, 0.3
    add_col_annotations(ax=ax, mdat=mdat, col_names=cols_to_plot, y0=y0, height=height, colormaps=colormaps)

    y0 = y0 + height + 0.01
    _, dendro = add_dendrogram_top(ax, Z, thres=thresh, y0=y0, height=0.6)

    # %% legend over all colormap categories (not just the plotted annotation columns)
    handles = []
    for col_name, colormap_dict in colormaps.items():
        for label, color in colormap_dict.items():
            handles.append(Line2D([0], [0], color=color, lw=4, label=f"{col_name}: {label}"))
    ax.legend(handles=handles, title="Sample annotation", bbox_to_anchor=(1.02, 1), loc="upper left")

    plot_path = save_dir / "figureS3a_stacked_barplot.pdf"
    ax.figure.savefig(plot_path, dpi=300, bbox_inches="tight")
    logger.info(f"saved Supplementary Figure 3a to {plot_path}")

    # %% dendrogram-leaf-color group assignment
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    mpl_colors = prop_cycle.by_key()["color"]
    cluster_colors = {f"C{i}": color for i, color in enumerate(mpl_colors)}
    cluster_colors["black"] = "#000000"

    df_tree = pd.DataFrame({"leaf_index": dendro["leaves"], "leaf_color_group": dendro["leaves_color_list"]})
    df_tree["leaf_color"] = df_tree["leaf_color_group"].map(cluster_colors)

    mdat = mdat.copy()
    mdat["leaf_color_group"] = df_tree["leaf_color_group"].values
    mdat["leaf_color"] = df_tree["leaf_color"].values

    props.reset_index().to_parquet(save_dir / "figureS3a_tma_composition.parquet")
    mdat.reset_index().to_parquet(save_dir / "figureS3a_metadata_with_dendrogram_colors.parquet")
    logger.info(f"saved Supplementary Figure 3a panels and data to {save_dir}")


if __name__ == "__main__":
    CLI(main)
