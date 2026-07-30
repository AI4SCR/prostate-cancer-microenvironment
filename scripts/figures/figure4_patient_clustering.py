# %%
"""Reproduce Figure 4a-c: patient-level cell-type composition clustering (P1-P6).

REPLACES an earlier Methods-text reconstruction of this panel (average-
linkage JSD clustering cut to a fixed k=6 via `scipy.fcluster(...,
criterion="maxclust")`) that was flagged as producing a confirmed mismatch
against the paper's reported (non-significant) KM result -- see
REPRODUCIBILITY.md's "CORRECTION" note on Figure 4a/c.

The actual generating script was since found: the old repo's
`000_paper/100_other_visualization/plot_stacked_frequencies.py`, run with
`var_name='label', group_var='pat_id'` (confirmed by the output filename it
produces, `metadata_with_dendrogram_colors_label_pat_id.parquet`, which
matches a file that's existed on shared storage all along -- see
missing_files.md). This is a 1:1 port of that script's `var_name='label'`
branch (paths only changed; see figure_script_mapping.md), and it is a
different, simpler algorithm than the Methods-text reconstruction:

1. Filter to cells from ROIs where `is_tumor == "yes"` only (not all ROIs).
2. Per-patient cell-type composition: pool ALL of a patient's cells
   directly (not per-core proportions max-pooled across cores) -- categorical
   value_counts with pseudocount 1, normalized, over all 35 labels
   (including "undefined" -- NOT excluded, unlike the earlier reconstruction).
3. Pairwise Jensen-Shannon distance via `scipy.spatial.distance.pdist(...,
   metric='jensenshannon')` -- natural-log JSD (scipy's default base), NOT
   base-2 like the earlier reconstruction used.
4. Average-linkage hierarchical clustering, cut by a HEIGHT threshold
   (`fcluster(Z, t=0.4, criterion="distance")`), not a fixed cluster count.

`compute_label_frequency`/`get_label_frequency_table` are ported inline from
the old repo's `datamodules/utils.py` (a small, self-contained utility, not
staged under LEGACY_DATA_DIR since it's two functions, not a whole module
worth importing).

Reads from `$EXPORT_DIR` (never `$BASE_DIR`). Writes the patient x cell-type
composition matrix, cluster assignments, and the clustered heatmap to
`$EXPORT_DIR/figures/figure4/`.
"""
from pathlib import Path

import pandas as pd
import seaborn as sns
from jsonargparse import CLI
from loguru import logger
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform

from prostate_cancer.utils import resolve_export_dir

DISTANCE_THRESHOLD = 0.4  # height cut, per plot_stacked_frequencies.py


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


def main(export_dir: Path | None = None):
    export_dir = export_dir or resolve_export_dir()
    save_dir = export_dir / "figures" / "figure4"
    save_dir.mkdir(parents=True, exist_ok=True)

    logger.info("loading exported tables")
    df_labels = pd.read_parquet(export_dir / "metadata.parquet").reset_index()
    clinical = pd.read_parquet(export_dir / "clinical.parquet")

    # %% FILTER TMAS TO INCLUDE: is_tumor == "yes" ROIs only
    df_sample_id = clinical.reset_index()
    df_sample_id = df_sample_id[~df_sample_id["is_tumor"].isna()]
    df_sample_id = df_sample_id[df_sample_id["is_tumor"] == "yes"]
    df_labels = df_labels.merge(df_sample_id, on="sample_id", how="inner")
    df_labels = df_labels.set_index(["sample_id", "object_id"])
    df_clusters = df_labels.copy()
    logger.info(f"{len(df_clusters)} cells from is_tumor=='yes' ROIs")

    # %% restrict clinical metadata to patients with at least one included ROI
    group_var = "pat_id"
    pat_cols = ["pat_id", "gs_grp", "os_status", "cause_of_death", "clinical_progr", "psa_progr", "disease_progr"]
    valid_tma_ids = df_clusters["tma_id"].unique().tolist()
    metadata = clinical[clinical["tma_id"].isin(valid_tma_ids)]
    metadata = metadata.set_index("tma_id")
    metadata = metadata[pat_cols]
    metadata = metadata.reset_index()
    metadata = metadata.drop_duplicates(subset=[group_var])
    metadata = metadata.set_index(group_var)

    # %% per-patient cell-type composition (pool all of a patient's cells directly)
    df_freqs = get_label_frequency_table(data=df_clusters, level="label", group_vars=[group_var])
    df_freqs, metadata = df_freqs.align(metadata, join="inner", axis=0)
    logger.info(f"patient-level composition matrix: {df_freqs.shape[0]} patients x {df_freqs.shape[1]} labels")

    # %% JSD hierarchical clustering, height cut (not a fixed k)
    dcond = pdist(df_freqs.values, metric="jensenshannon")
    Z = linkage(dcond, method="average")
    fixed_clusters = fcluster(Z, t=DISTANCE_THRESHOLD, criterion="distance")
    patient_clusters = pd.Series([f"P{c}" for c in fixed_clusters], index=df_freqs.index, name="patient_cluster")
    logger.info(f"cluster sizes (distance threshold {DISTANCE_THRESHOLD}):\n{patient_clusters.value_counts().sort_index()}")

    composition = df_freqs
    composition.reset_index().to_parquet(save_dir / "figure4a_patient_composition.parquet")
    patient_clusters.reset_index().to_parquet(save_dir / "figure4a_patient_clusters.parquet")

    # %% Figure 4a: clustered heatmap, rows ordered/colored by patient cluster
    n_clusters = patient_clusters.nunique()
    cluster_colors = dict(zip(sorted(patient_clusters.unique()), sns.color_palette("tab10", n_clusters)))
    row_colors = patient_clusters.map(cluster_colors)
    cg = sns.clustermap(
        composition,
        row_linkage=Z,
        col_cluster=True,
        row_colors=row_colors,
        figsize=(14, 16),
        cmap="viridis",
        yticklabels=False,
    )
    cg.ax_heatmap.set_title(f"Figure 4a -- patient-level cell-type composition ({n_clusters} clusters)")
    cg.figure.savefig(save_dir / "figure4a_composition_heatmap.png", dpi=200)

    logger.info(f"saved figure 4a panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
