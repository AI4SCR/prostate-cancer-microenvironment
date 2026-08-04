# %%
"""Reproduce Figure 4b: mean cellular metagroup distribution per patient cluster.

1:1 port of the newly-pulled old repo's
`000_paper/sync_paper/05-heterogeneity/stacked-frequencies-label.py`
(`group_var='pat_id'` branch, its second half -- lines ~322-367, the block
gated `if group_var == 'pat_id':`). See figure_script_mapping.md.

Averages per-patient cell-type (`label`) composition within each patient
cluster, maps each of the 34 fine-grained labels to one of 17 coarser
"metagroup" lineages via `resources/metalabels.yaml`'s `metagroups` key
(ported from `000_paper/metalabels.yaml`), then sums within each
(cluster, metagroup) pair to get the average metagroup composition per
cluster.

Two disclosed deviations, both consistent with the already-established
convention in figure4_survival.R (see that script's docstring for the full
rationale, not repeated here):
- Reads `figure4_patient_clustering.py`'s own composition matrix
  (`figure4a_patient_composition.parquet`) rather than recomputing
  per-patient label frequencies inline -- identical computation
  (`get_label_frequency_table`), already verified byte-equivalent between
  the two legacy script versions, just not duplicated a third time.
- Uses the precomputed patient-cluster partition
  ($LEGACY_DATA_DIR/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet,
  `leaf_color_group`, excluding "black") instead of live-recomputing the
  dendrogram/leaf-coloring here, for the same reason figure4_survival.R
  does: our own fcluster reconstruction's near-singleton leaves don't
  exactly match matplotlib's dendrogram(color_threshold=...) leaf-coloring
  that produced this file.

Reads from `$EXPORT_DIR`/`$LEGACY_DATA_DIR`/`resources/`. Writes to
`$OUTPUT_FIGURES_DIR/figure4/`.
"""
from pathlib import Path

import matplotlib
import pandas as pd
import yaml
from jsonargparse import CLI
from loguru import logger
from matplotlib import pyplot as plt

from prostate_cancer.utils import resolve_export_dir, resolve_legacy_dir, resolve_output_figures_dir

matplotlib.use("Agg")


def main(export_dir: Path | None = None, legacy_dir: Path | None = None):
    export_dir = export_dir or resolve_export_dir()
    legacy_dir = legacy_dir or resolve_legacy_dir()
    save_dir = resolve_output_figures_dir() / "figure4"
    save_dir.mkdir(parents=True, exist_ok=True)
    resources_dir = Path(__file__).resolve().parents[2] / "resources"

    df_freqs = pd.read_parquet(save_dir / "figure4a_patient_composition.parquet").set_index("pat_id")

    patient_clusters = pd.read_parquet(
        legacy_dir / "5-niches" / "barplot_data" / "metadata_with_dendrogram_colors_label_pat_id.parquet"
    )
    patient_clusters = patient_clusters[patient_clusters["leaf_color_group"] != "black"]

    df_freqs, patient_clusters = df_freqs.align(patient_clusters, join="inner", axis=0)
    df_freqs = df_freqs.copy()
    df_freqs["grouping"] = patient_clusters["leaf_color_group"].values
    logger.info(f"{len(df_freqs)} patients across {df_freqs['grouping'].nunique()} clusters")

    # %% average per-label composition within each cluster
    df_averaged = df_freqs.groupby("grouping").mean()
    df_averaged_long = df_averaged.reset_index().melt(id_vars="grouping", var_name="label", value_name="average_proportion")

    # %% map labels to metagroups, sum within (cluster, metagroup)
    with open(resources_dir / "metalabels.yaml") as f:
        lineage_map = yaml.safe_load(f)["metagroups"]
    df_averaged_long["metacelltype"] = df_averaged_long["label"].map(lineage_map)
    df_averaged_long_meta = df_averaged_long.groupby(["grouping", "metacelltype"]).sum(numeric_only=True).reset_index()
    df_freqs_meta_wide = df_averaged_long_meta.pivot(index="grouping", columns="metacelltype", values="average_proportion")

    # %% plot
    with open(resources_dir / "colormaps.yaml") as f:
        colormaps = yaml.safe_load(f)
    colordict_meta = colormaps["metagroups"]

    prop_cycle = plt.rcParams["axes.prop_cycle"]
    mpl_colors = prop_cycle.by_key()["color"]
    cluster_colors = {f"C{i}": color for i, color in enumerate(mpl_colors)}
    cluster_colors["black"] = "#000000"

    fig, ax = plt.subplots(figsize=(14, 12))
    df_freqs_meta_wide.plot(
        kind="bar",
        stacked=True,
        color=[colordict_meta[c] for c in df_freqs_meta_wide.columns],
        ax=ax,
    )
    # "C1".."C6" -> "P1".."P6" for the displayed labels only: the published
    # figure labels clusters "P1"-"P6", not the precomputed file's own
    # "leaf_color_group" values. Underlying grouping/color-key values are
    # left as "C1".."C6" (matplotlib's own C0-C9 cycle convention, which the
    # cluster_colors dict above is keyed by) -- only the tick text changes.
    # No script in the legacy repo performs this relabeling; applied here to
    # match the published figure (same fix as figure4c_patient_cluster_km.R).
    display_labels = [g.replace("C", "P", 1) for g in df_freqs_meta_wide.index]
    ax.set_xticklabels(display_labels, rotation=0)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Average Proportion")
    ax.set_title("Figure 4b -- average cell type proportions by cluster")
    for tick, group in zip(ax.get_xticks(), df_freqs_meta_wide.index):
        ax.add_patch(
            plt.Rectangle(
                (tick - 0.4, -0.05), 0.8, 0.05, color=cluster_colors[group], transform=ax.get_xaxis_transform(), clip_on=False
            )
        )
    ax.legend(title="Cell Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    save_path = save_dir / "figure4b_metagroup_barplot.pdf"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    logger.info(f"saved Figure 4b to {save_path}")


if __name__ == "__main__":
    CLI(main)
