"""Reproduce Supplementary Figure 2a-c: compartment-specific UMAPs (immune/
epithelial/endothelial), colored by cell type and marker intensity.

1:1 port of `archive/scripts/02-umaps/0-umaps-main-types.py`'s plotting
loops. Unlike `figure2_umap.py` (which subsets the single all-cells UMAP
embedding by `main_group` for its own separate per-compartment panels),
the legacy pipeline fit a genuinely SEPARATE UMAP per compartment
(`main_group=immune`/`epithelial`/`endothelial`, its own `reducer.pkl`
each) -- confirmed via `archive/scripts/02-umaps/0-umaps-main-types.py`'s
`params = list(product(..., main_groups, ...))` and per-`main_group`
`reducer_path`. `figure2_umap.py`'s compartment-subset panels are a
different thing (Figure 2b's own main_group breakdown, not this
supplementary figure) and are not touched by this script.

Each compartment's embedding is ported from its own legacy `reducer.pkl`
via `scripts/port/port_umap_reducer.py` (that script's own pinned pixi env,
NOT this repo's `.venv` -- see its docstring):

    cd scripts/port
    for mg in immune epithelial endothelial; do
      pixi run python port_umap_reducer.py \\
        "/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/2-umaps/2-main_groups/n_neighbors=50-min_dist=0.1-engine=umap-learn-main_group=${mg}-excl_markers=dna1_dna2_fap_icsk1_icsk2_icsk3/reducer.pkl" \\
        "../../data/figures/figureS2_main_groups_umap/${mg}/umap_embeddings.parquet"
    done

(legacy also fit `stromal`/`undefined` compartments, not part of this
supplementary figure per its legend, so not ported here.)
"""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from ai4bmr_learn.utils.sampling import sample_min_per_group_then_uniform
from jsonargparse import CLI
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_rgba

from prostate_cancer.utils import get_colormap_dict, resolve_export_dir, resolve_output_figures_dir

NUM_SAMPLES_PER_PLOT = 100_000
MAIN_GROUPS = ["immune", "epithelial", "endothelial"]


def plot_points(
    *, reducer_embedding, subset_points=None, values=None, cmap=None, cbar: bool = False,
    labels=None, color_key=None, shuffle: bool = True,
):
    fig, ax = plt.subplots()
    x = reducer_embedding[:, 0]
    y = reducer_embedding[:, 1]

    if subset_points is not None:
        x = x[subset_points]
        y = y[subset_points]

    if labels is not None:
        l = labels[subset_points]
        c = np.array([to_rgba(color_key[i], alpha=1) for i in l])
    elif values is not None:
        v = values[subset_points]
        c = cmap(v)
    else:
        raise ValueError("Must provide labels or values")

    if shuffle:
        order = np.arange(len(x))
        np.random.shuffle(order)
        x, y = x[order], y[order]
        c = c[order]

    ax.scatter(x=x, y=y, s=1, c=c)
    ax.set_axis_off()

    if values is not None and cmap is not None and cbar:
        sm = ScalarMappable(cmap=cmap)
        sm.set_array(values[subset_points])
        fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.046, pad=0.04)

    return ax


def main(export_dir: Path | None = None):
    export_dir = export_dir or resolve_export_dir()
    save_dir = resolve_output_figures_dir() / "figureS2"
    save_dir.mkdir(parents=True, exist_ok=True)

    logger.info("loading exported tables")
    metadata = pd.read_parquet(export_dir / "metadata.parquet")
    data = pd.read_parquet(export_dir / "intensity_normalized.parquet")
    metadata, data = metadata.align(data, axis=0, join="inner")

    ported_dir = Path(__file__).resolve().parents[2] / "data" / "figures" / "figureS2_main_groups_umap"

    for main_group in MAIN_GROUPS:
        reducer_path = ported_dir / main_group / "umap_embeddings.parquet"
        assert reducer_path.exists(), (
            f"{reducer_path} missing -- run scripts/port/port_umap_reducer.py "
            f"for main_group={main_group} first, see this script's docstring"
        )
        embedding_df = pd.read_parquet(reducer_path).set_index(["sample_id", "object_id"])
        index = embedding_df.index
        embedding = embedding_df[["umap_1", "umap_2"]].values

        group_dir = save_dir / main_group
        group_dir.mkdir(parents=True, exist_ok=True)

        # %% label panel
        logger.info(f"[{main_group}] plotting UMAP colored by label")
        colormap_dict = get_colormap_dict(name="label")
        md = metadata.loc[index].copy()
        labels = md["label"].values

        grouped = md.groupby("label", observed=True)
        min_per_group = min(NUM_SAMPLES_PER_PLOT // grouped.ngroups, grouped.size().min())
        md_sub = sample_min_per_group_then_uniform(grouped=grouped, n=NUM_SAMPLES_PER_PLOT, min_per_group=min_per_group, random_state=0)

        subset_points = np.zeros(len(index), dtype=bool)
        subset_points[index.isin(md_sub.index)] = True

        ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, labels=labels, color_key=colormap_dict)
        ax.figure.tight_layout()
        ax.figure.savefig(group_dir / "label=label.pdf", transparent=True)
        plt.close(ax.figure)

        # %% per-marker intensity panels
        for value in data.columns:
            logger.info(f"[{main_group}] plotting UMAP colored by intensity of {value}")
            values = data.loc[index, value].values
            cmap = matplotlib.colormaps["Reds"]
            norm = Normalize(vmin=values.min(), vmax=values.max())

            grouped = md.groupby("label", observed=True)
            min_per_group = min(NUM_SAMPLES_PER_PLOT // grouped.ngroups, grouped.size().min())
            md_sub = sample_min_per_group_then_uniform(grouped=grouped, n=NUM_SAMPLES_PER_PLOT, min_per_group=min_per_group, random_state=0)

            subset_points = np.zeros(len(index), dtype=bool)
            subset_points[index.isin(md_sub.index)] = True

            ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, values=norm(values), cmap=cmap, cbar=True, shuffle=False)
            ax.figure.tight_layout()
            ax.figure.savefig(group_dir / f"value={value}.pdf", transparent=True)
            plt.close(ax.figure)

        logger.info(f"[{main_group}] saved all panels to {group_dir}")

    logger.info(f"saved all Supplementary Figure 2a-c panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
