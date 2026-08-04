# %%
"""Reproduce Figure 6c's "top" panel: mean cell-type composition across
cores of niches 6, 16, 17, and 18.

1:1 port of the old repo's
`000_paper/sync_paper/06-spatial-niches/composition/niche_composition.py`'s
"MEAN COMPOSITION BARPLOT" section, same as `figureS4b_niche_mean_composition.py`
(see that script's docstring for the shared source and the parts of
`niche_composition.py` NOT ported), but with the `niche_order` list
restricted to the 4 niches this panel names -- a disclosed, parameter-only
filter of the same verbatim code, not a logic change (same permitted
category as e.g. `figure6_km_niche6.R` scanning a different niche each
time).

Only the "top" mean-composition panel is produced here. Figure 6c's
"bottom" panel ("two example cores") is NOT code-derived -- no legacy
script plots per-core composition bars, only the corpus-wide mean shown
here -- and is flagged as a gap in open-questions.md rather than attempted.

Reads `clusters_annotated_v2.parquet` from `LEGACY_DATA_DIR` and
`resources/colormaps.yaml`. Writes to `$OUTPUT_FIGURES_DIR/figure6/`.
"""
from pathlib import Path

import matplotlib
import pandas as pd
import yaml
from jsonargparse import CLI
from loguru import logger
from matplotlib import pyplot as plt

from prostate_cancer.utils import resolve_legacy_dir, resolve_output_figures_dir

matplotlib.use("Agg")

# Restricted to niches 6/16/17/18 (Figure 6c), per the same fixed niche_order
# list used throughout the legacy niche-visualization scripts -- see
# figureS4b_niche_mean_composition.py's NICHE_ORDER for the full 18-niche list.
NICHE_ORDER = [
    "tumorERG+p53+_ProlifLuminal",  # niche 6
    "immune_bloodvessels_CAF1(CD105-)",  # niche 16
    "Macrophages_Tcells_CAF1(CD105-)",  # niche 17
    "TLS",  # niche 18
]
MIN_CELLS_PER_NICHE_SAMPLE = 15


def main(legacy_dir: Path | None = None):
    legacy_dir = legacy_dir or resolve_legacy_dir()
    save_dir = resolve_output_figures_dir() / "figure6"
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
    df_comp_mean = df_composition.reset_index()
    df_comp_mean = df_comp_mean.drop(columns=[group_var]).groupby(niche_col).mean()

    with open(resources_dir / "colormaps.yaml") as f:
        color_maps = yaml.safe_load(f)
    color_dict_label = color_maps["label"]

    df_comp_mean = df_comp_mean.reindex(NICHE_ORDER)
    colors = [color_dict_label[c] for c in df_comp_mean.columns]

    ax = df_comp_mean.plot(kind="bar", stacked=True, figsize=(10, 6), color=colors)

    plt.title("Mean Cell Type Composition per Niche (niches 6, 16, 17, 18)")
    plt.ylabel("Mean Composition")
    plt.xlabel("Niche")
    plt.legend(title="Cell Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    fig_path = save_dir / "figure6c_niche_composition_barplot.pdf"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    logger.info(f"saved Figure 6c top panel to {fig_path}")


if __name__ == "__main__":
    CLI(main)
