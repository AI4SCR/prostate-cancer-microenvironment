# %%
"""Reproduce Figure 2b: UMAP of all cells colored by compartment/cell type/patient/markers.

1:1 port of the old repo's `000_paper/02_umaps/0-umaps.py` (see
`figure_script_mapping.md`). That script explores a `params` sweep
(`engine` x `n_neighbors` x `min_dist` x `exclude_markers`) but only ever
actually runs ONE active configuration (line 27 of the original --
everything else is commented out): `engine='umap-learn'`, `n_neighbors=50`,
`min_dist=0.1`, excluding the same non-marker channels this repo calls
`NON_MARKER_CHANNELS`. That config is hardcoded here rather than ported as
a sweep, since the sweep was never actually exercised.

Differences from the legacy script, all disclosed:
- Loads `intensity_normalized.parquet` from `$EXPORT_DIR` directly, instead
  of loading raw intensities via a live `ai4bmr_datasets.PCa()` call and
  then calling `normalize(data, exclude_zeros=True)`. This repo's
  `intensity_normalized.parquet` (via `export_for_r.py`) IS
  `normalize(raw_intensity, exclude_zeros=True)` -- confirmed byte-identical
  to the legacy export (see REPRODUCIBILITY.md) -- so this is the same
  values, not a different computation, just reading a cache of it that
  already exists instead of recomputing it.
- Legacy computes the UMAP fit on the FULL dataset (no subsampling before
  `reducer.fit()`) and only subsamples for the scatter plots afterward, via
  `sample_min_per_group_then_uniform`. Preserved here exactly -- this is a
  slow, full-dataset UMAP.fit() on ~2.19M cells, not a quick 50k-cell demo.

Reads from `$EXPORT_DIR` (never `$BASE_DIR`). Writes the fitted embedding
and all plot panels to `$EXPORT_DIR/figures/figure2/`.
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

from prostate_cancer.utils import NON_MARKER_CHANNELS, get_colormap_dict, resolve_export_dir

N_NEIGHBORS = 50
MIN_DIST = 0.1
NUM_SAMPLES_PER_PLOT = 100_000  # legacy's num_samples for sample_min_per_group_then_uniform


def plot_points(
    *, reducer_embedding, subset_points=None, values=None, cmap=None, cbar: bool = False,
    labels=None, color_key=None, all_points: bool = False, shuffle: bool = True, order=None,
):
    fig, ax = plt.subplots()
    if all_points:
        x = reducer_embedding[:, 0]
        y = reducer_embedding[:, 1]
        c = to_rgba("#D3D3D3", alpha=0.1)
        ax.scatter(x=x, y=y, s=0.1, color=c)

    if order is not None:
        x = reducer_embedding[order, 0]
        y = reducer_embedding[order, 1]
    else:
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


def main(export_dir: Path | None = None, seed: int = 0):
    import umap

    export_dir = export_dir or resolve_export_dir()
    save_dir = export_dir / "figures" / "figure2"
    save_dir.mkdir(parents=True, exist_ok=True)

    # %% load exported cells (== ds.intensity/ds.metadata + normalize(exclude_zeros=True), see docstring)
    logger.info("loading exported tables")
    metadata = pd.read_parquet(export_dir / "metadata.parquet")
    data = pd.read_parquet(export_dir / "intensity_normalized.parquet")
    metadata, data = metadata.align(data, axis=0, join="inner")
    clinical = pd.read_parquet(export_dir / "clinical.parquet")
    sid_to_pid = clinical["pat_id"].to_dict()

    exclude_markers = NON_MARKER_CHANNELS
    fit_data = data.loc[:, ~data.columns.isin(exclude_markers)].copy()

    # %% COMPUTE UMAP on the full dataset (no subsampling for the fit itself)
    reducer_path = save_dir / "reducer_embedding.parquet"
    if reducer_path.exists():
        logger.info(f"reusing existing embedding at {reducer_path}")
        embedding_df = pd.read_parquet(reducer_path)
        index = embedding_df.index
        embedding = embedding_df[["umap_1", "umap_2"]].values
    else:
        logger.info(f"computing UMAP for {len(fit_data)} cells, n_neighbors={N_NEIGHBORS}, min_dist={MIN_DIST}, excluding {exclude_markers}")
        reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, metric="euclidean", random_state=seed)
        reducer.fit(fit_data.values)
        index = fit_data.index
        embedding = reducer.embedding_
        embedding_df = pd.DataFrame(embedding, index=index, columns=["umap_1", "umap_2"])
        embedding_df.to_parquet(reducer_path)
        logger.info(f"saved embedding to {reducer_path}")

    # %% label/main_group/pat_id panels
    for label in ["main_group", "label", "pat_id"]:
        logger.info(f"plotting UMAP colored by {label}")
        colormap_dict = get_colormap_dict(name=label)

        md = metadata.loc[index].copy()
        md["pat_id"] = md.index.get_level_values("sample_id").map(sid_to_pid)
        labels = md[label].values

        grouped = md.groupby(label, observed=True)
        min_per_group = min(NUM_SAMPLES_PER_PLOT // grouped.ngroups, grouped.size().min())
        md_sub = sample_min_per_group_then_uniform(grouped=grouped, n=NUM_SAMPLES_PER_PLOT, min_per_group=min_per_group, random_state=0)

        subset_points = np.zeros(len(index), dtype=bool)
        subset_points[index.isin(md_sub.index)] = True

        ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, labels=labels, color_key=colormap_dict)
        ax.figure.tight_layout()
        ax.figure.savefig(save_dir / f"label={label}.pdf", transparent=True)
        plt.close(ax.figure)

    # %% per-main_group subset panels
    for main_group in ["immune", "epithelial", "endothelial", "stromal"]:
        logger.info(f"plotting UMAP subset for main_group={main_group}")
        colormap_dict = get_colormap_dict(name="label")

        md = metadata.loc[index].copy()
        labels = md["label"].values

        md_group = md[md["main_group"] == main_group]
        grouped = md_group.groupby("label", observed=True)
        min_per_group = min(NUM_SAMPLES_PER_PLOT // grouped.ngroups, grouped.size().min())
        md_sub = sample_min_per_group_then_uniform(grouped=grouped, n=NUM_SAMPLES_PER_PLOT, min_per_group=min_per_group, random_state=0)

        subset_points = np.zeros(len(index), dtype=bool)
        subset_points[index.isin(md_sub.index)] = True

        ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, labels=labels, color_key=colormap_dict, all_points=True)
        ax.figure.tight_layout()
        ax.figure.savefig(save_dir / f"label={main_group}.pdf", transparent=True)
        plt.close(ax.figure)

    # %% per-marker intensity panels
    for value in data.columns:
        logger.info(f"plotting UMAP colored by intensity of {value}")
        md = metadata.loc[index].copy()
        md["pat_id"] = md.index.get_level_values("sample_id").map(sid_to_pid)

        values = data.loc[index, value].values
        cmap = matplotlib.colormaps["Reds"]
        norm = Normalize(vmin=values.min(), vmax=values.max())

        grouped = md.groupby("label", observed=True)
        min_per_group = min(NUM_SAMPLES_PER_PLOT // grouped.ngroups, grouped.size().min())
        md_sub = sample_min_per_group_then_uniform(grouped=grouped, n=NUM_SAMPLES_PER_PLOT, min_per_group=min_per_group, random_state=0)

        subset_points = np.zeros(len(index), dtype=bool)
        subset_points[index.isin(md_sub.index)] = True

        ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, values=norm(values), cmap=cmap, cbar=True, shuffle=False, order=None)
        ax.figure.tight_layout()
        ax.figure.savefig(save_dir / f"value={value}.pdf", transparent=True)
        plt.close(ax.figure)

    logger.info(f"saved all Figure 2b panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
