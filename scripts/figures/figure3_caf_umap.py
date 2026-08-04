# %%
"""Reproduce Figure 3a: UMAP of CAF cells colored by subcluster.

1:1 port of the old repo's `000_paper/02_umaps/0-umaps-cafs.py` (see
`figure_script_mapping.md`). Legacy's active `params` sweep has TWO
configs, not one (both n_neighbors=50, min_dist=0.1, engine='umap-learn'):
(a) excluding only the standard non-marker channels (`NON_MARKER_CHANNELS`)
-- i.e. essentially all markers, and (b) excluding every marker except the
9 CAF-relevant ones the paper's Methods names. Both are computed here;
config (b) is the one this repo's earlier notes identified as producing the
actual Fig 3a panel ("UMAP ... computed using CAF-relevant markers").

Cell filter: legacy selects `metadata.label.str.contains('CAF')` -- NOT
`main_group == 'stromal'` (stromal includes non-CAF cells like pericytes;
an earlier, non-faithful version of this script used the broader filter --
see discrepancies.md).

Normalization: legacy calls `normalize(df, exclude_zeros=True)` on the
CAF-filtered, marker-subsetted RAW data, per config -- i.e. the min-max/
censoring statistics are computed on that specific filtered population and
marker subset, not globally. This is NOT the same computation as this
repo's `intensity_normalized.parquet` (which normalizes ALL cells x ALL
markers together), so that cache can't be reused here -- this script loads
`intensity.parquet` (raw) instead and normalizes it itself, matching legacy.

Both UMAP fits (data-loading/export mechanics only, not logic) now live in
`scripts/data/figure3_caf_umap_embedding.py`. This script only reads those
cached embeddings and plots -- it still loads and normalizes the raw
CAF-filtered intensity table itself, since the intensity panels need it for
coloring (a cheap, non-UMAP step) independent of the cached embeddings.

Reads from `$EXPORT_DIR`, requires
`scripts/data/figure3_caf_umap_embedding.py` to have already produced
`$EXPORT_DIR/figures/figure3/{config_name}/umap_embedding.parquet` for both
configs. Writes all plot panels to `$EXPORT_DIR/figures/figure3/`.
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

from prostate_cancer.utils import NON_MARKER_CHANNELS, get_colormap_dict, normalize, resolve_export_dir, resolve_output_figures_dir

ALL_MARKERS = [
    "smooth_muscle_actin", "prostate_specific_antigen", "vimentin", "collagen1", "synaptophysin", "keratin5",
    "yap1", "pan_keratin", "ces1", "egr1", "cd31", "cd45", "cd44", "fap", "fox_p3", "cd4", "e_cadherin", "cd68",
    "cd66b", "cd20", "cd8a", "cd11b", "p63", "beta_catenin", "pdpn", "cd105", "ki_67", "p53", "cd3", "erg",
    "c_casp3", "cnn1", "keratin8_18", "cd146", "ar", "dna1", "dna2", "icsk1", "icsk2", "icsk3",
]
CAF_MARKERS = {"smooth_muscle_actin", "vimentin", "collagen1", "cd146", "cnn1", "cd105", "ar", "egr1", "ces1"}
CONFIGS = {
    "excl_markers": sorted(NON_MARKER_CHANNELS),
    "caf_markers_only": sorted(set(ALL_MARKERS) - CAF_MARKERS),
}
NUM_SAMPLES_PER_PLOT = 100_000


def plot_points(
    *, reducer_embedding, subset_points=None, values=None, cmap=None, cbar: bool = False,
    labels=None, color_key=None, all_points: bool = False, shuffle: bool = True,
):
    fig, ax = plt.subplots()
    if all_points:
        x = reducer_embedding[:, 0]
        y = reducer_embedding[:, 1]
        c = to_rgba("#D3D3D3", alpha=0.1)
        ax.scatter(x=x, y=y, s=0.1, color=c)

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
    save_dir = resolve_output_figures_dir() / "figure3"
    save_dir.mkdir(parents=True, exist_ok=True)

    logger.info("loading exported tables")
    metadata = pd.read_parquet(export_dir / "metadata.parquet")
    raw_intensity = pd.read_parquet(export_dir / "intensity.parquet")
    metadata, raw_intensity = metadata.align(raw_intensity, axis=0, join="inner")
    clinical = pd.read_parquet(export_dir / "clinical.parquet")
    sid_to_pid = clinical["pat_id"].to_dict()

    caf_filter = metadata["label"].str.contains("CAF")
    data = raw_intensity.loc[caf_filter, :]
    caf_metadata = metadata.loc[caf_filter, :]
    logger.info(f"{len(data)} CAF cells (label contains 'CAF')")

    for config_name in CONFIGS:
        config_dir = save_dir / config_name
        config_dir.mkdir(parents=True, exist_ok=True)

        embedding_path = config_dir / "umap_embedding.parquet"
        assert embedding_path.exists(), f"{embedding_path} missing -- run scripts/data/figure3_caf_umap_embedding.py first"
        embedding_df = pd.read_parquet(embedding_path)
        index = embedding_df.index
        embedding = embedding_df[["umap_1", "umap_2"]].values

        # %% label/main_group/pat_id panels
        for label in ["main_group", "label", "pat_id"]:
            colormap_dict = get_colormap_dict(name=label)
            md = caf_metadata.loc[index].copy()
            md["pat_id"] = md.index.get_level_values("sample_id").map(sid_to_pid)
            labels = md[label].values

            grouped = md.groupby(label, observed=True)
            min_per_group = min(NUM_SAMPLES_PER_PLOT // grouped.ngroups, grouped.size().min())
            md_sub = sample_min_per_group_then_uniform(grouped=grouped, n=NUM_SAMPLES_PER_PLOT, min_per_group=min_per_group, random_state=0)

            subset_points = np.zeros(len(index), dtype=bool)
            subset_points[index.isin(md_sub.index)] = True

            ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, labels=labels, color_key=colormap_dict)
            ax.figure.tight_layout()
            ax.figure.savefig(config_dir / f"label={label}.pdf", transparent=True)
            plt.close(ax.figure)

        # %% intensity panels
        df_norm = normalize(data.loc[index, :], exclude_zeros=True)
        for value in df_norm.columns:
            md = caf_metadata.loc[index].copy()
            md["pat_id"] = md.index.get_level_values("sample_id").map(sid_to_pid)
            values = df_norm[value].values
            cmap = matplotlib.colormaps["Reds"]
            norm = Normalize(vmin=values.min(), vmax=values.max())

            grouped = md.groupby("label", observed=True)
            min_per_group = min(NUM_SAMPLES_PER_PLOT // grouped.ngroups, grouped.size().min())
            md_sub = sample_min_per_group_then_uniform(grouped=grouped, n=NUM_SAMPLES_PER_PLOT, min_per_group=min_per_group, random_state=0)

            subset_points = np.zeros(len(index), dtype=bool)
            subset_points[index.isin(md_sub.index)] = True

            ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, values=norm(values), cmap=cmap, cbar=True, shuffle=False)
            ax.figure.tight_layout()
            ax.figure.savefig(config_dir / f"value={value}.pdf", transparent=True)
            plt.close(ax.figure)

        # %% per-CAF-subtype subset panels
        colormap_dict = get_colormap_dict(name="label")
        for cell_type in sorted(set(filter(lambda x: "CAF" in x, metadata["label"].unique()))):
            md = caf_metadata.loc[index].copy()
            labels = md["label"].values
            select = md["label"] == cell_type

            subset_points = np.zeros(len(index), dtype=bool)
            subset_points[index.isin(md.index[select])] = True

            ax = plot_points(reducer_embedding=embedding, subset_points=subset_points, labels=labels, color_key=colormap_dict)
            ax.figure.tight_layout()
            ax.figure.savefig(config_dir / f"cell_type={cell_type}.pdf", transparent=True)
            plt.close(ax.figure)

        logger.info(f"[{config_name}] saved all panels to {config_dir}")

    logger.info(f"saved all Figure 3a panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
