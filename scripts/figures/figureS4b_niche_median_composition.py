"""Reproduce Supplementary Figure 4b: median cell type composition per niche.

Port of the "CORRECT VERSION WITH MEDIAN INSTEAD OF MEAN" block in
`sync_paper/06-spatial-niches/composition/niche_composition.py` (as of
legacy commit `84c1f2d`, "adapted niche composition figure construction" --
newer than the commit `figureS4b_niche_mean_composition.py` was ported
from, which only computed `df_comp_median` without ever plotting it). The
paper's published S4b legend describes median-per-cell-type-across-samples
composition, rescaled per niche to sum to 1 -- not the mean composition
`figureS4b_niche_mean_composition.py` produces; that script is left
unmodified (see `figures.md`'s S4b Issue note) and this is a separate,
dedicated script for the median variant.
"""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import yaml
from jsonargparse import CLI
from loguru import logger
from matplotlib import pyplot as plt

from prostate_cancer.utils import resolve_legacy_dir, resolve_output_figures_dir

matplotlib.use("Agg")

NICHE_ORDER = [
    "luminal",
    "luminal_infiltrated",
    "basal_luminal_glands",
    "tumor(ERG+)",
    "tumor(ERG+)_luminal",
    "tumorERG+p53+_ProlifLuminal",
    "tumor_CAF1(CD105High)",
    "tumor_CAF1_lymphocytes",
    "luminal_CAF1(CD105High)",
    "luminal_CAF1(CD105-)",
    "bloodvessels_CAF1(CD105-)",
    "CAF1_CD105-_infiltrated",
    "CAFs_lymphocytes",
    "CAF2s_enriched",
    "CAF2(AR-)_enriched",
    "immune_bloodvessels_CAF1(CD105-)",
    "Macrophages_Tcells_CAF1(CD105-)",
    "TLS",
]
MIN_CELLS_PER_NICHE_SAMPLE = 15


def main(legacy_dir: Path | None = None):
    legacy_dir = legacy_dir or resolve_legacy_dir()
    save_dir = resolve_output_figures_dir() / "figureS4"
    save_dir.mkdir(parents=True, exist_ok=True)
    resources_dir = Path(__file__).resolve().parents[2] / "resources"

    celltype_col = "label"
    niche_col = "niche"
    group_var = "sample_id"
    cluster_path = legacy_dir / "5-niches" / "annotation" / "clusters_annotated_v2.parquet"
    df_clusters = pd.read_parquet(cluster_path)
    logger.info(f"clusters shape: {df_clusters.shape}")

    df_counts_per_niche = df_clusters.groupby([niche_col, group_var])[celltype_col].value_counts().rename("count").unstack(fill_value=0)
    df_niche_size = df_clusters.groupby([niche_col, group_var]).size().rename("niche_size")
    # don't consider niches with very few cells
    df_niche_size[df_niche_size <= MIN_CELLS_PER_NICHE_SAMPLE] = 0
    assert df_counts_per_niche.index.equals(df_niche_size.index), "Indices do not match!"
    df_composition = df_counts_per_niche.div(df_niche_size, axis=0).fillna(0)
    df_composition = df_composition[df_niche_size > 0]
    df_comp_median = df_composition.reset_index()
    df_comp_median = df_comp_median.drop(columns=[group_var]).groupby(niche_col).median()

    with open(resources_dir / "colormaps.yaml") as f:
        color_maps = yaml.safe_load(f)
    color_dict_label = color_maps["label"]

    df_comp_median = df_comp_median.reindex(NICHE_ORDER)

    # Rescale each niche composition so it sums to 1
    row_sums = df_comp_median.sum(axis=1)
    df_comp_median = df_comp_median.div(row_sums.replace(0, np.nan), axis=0).fillna(0)

    colors = [color_dict_label[c] for c in df_comp_median.columns]

    ax = df_comp_median.plot(
        kind="bar",
        stacked=True,
        figsize=(10, 6),
        color=colors,
    )

    plt.title("Median Cell Type Composition per Niche")
    plt.ylabel("Rescaled Median Composition")
    plt.xlabel("Niche")
    plt.ylim(0, 1)
    plt.legend(
        title="Cell Type",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    plt.tight_layout()
    fig_path = save_dir / "figureS4b_niche_composition_barplot_median.pdf"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    logger.info(f"saved Supplementary Figure 4b (median variant) to {fig_path}")


if __name__ == "__main__":
    CLI(main)
