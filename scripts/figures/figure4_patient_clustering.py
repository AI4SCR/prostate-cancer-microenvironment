# %%
"""Reproduce Figure 4a-b: patient-level cell-type composition clustering.

No prior script in this repo computed this -- reconstructed from the paper's
Methods ("Cell Type Proportion Quantification"):

1. Cell-type proportions per physical core (`tma_id`, not `sample_id` --
   interrupted-acquisition duplicates must be pooled first, see
   REPRODUCIBILITY.md), over the 34 annotated cell types (excluding
   "undefined"), with a pseudocount of 1 added to every type before
   normalizing (avoids zeros for the log-ratio/JSD steps).
2. Patient-level composition vector via max-pooling: for each cell type, the
   maximum proportion observed across that patient's cores.
3. Hierarchical clustering (average linkage) of patient compositions using
   Jensen-Shannon divergence as the distance metric, cut to 6 clusters
   (P1-P6, matching the paper's reported "six patient groups").

Reads from `$EXPORT_DIR` (never `$BASE_DIR` -- see REPRODUCIBILITY.md).
Writes the patient x cell-type composition matrix, cluster assignments, and
the clustered heatmap to `$EXPORT_DIR/figures/figure4/`.

CONFIRMED MISMATCH, unresolved: against the real dataset this produces P1-P6
cluster sizes 3/101/23/2/65/1 -- three near-singleton clusters. The paper
reports the resulting KM analysis (Fig 4c) as NOT significant; this
reconstruction's clusters give log-rank p=6.65e-09, driven by those tiny
clusters' volatile survival curves. The Methods text doesn't specify how the
dendrogram was cut into exactly 6 groups (a fixed height threshold, not used
here, would plausibly give a more balanced split than forcing k=6). See
REPRODUCIBILITY.md Known discrepancies. Figure 4d-e (Cox PH on the same
underlying composition, independent of this clustering step) DO reproduce
the paper's result exactly, so the composition/CLR computation itself is
validated -- only the P1-P6 cut point is in question.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from jsonargparse import CLI
from loguru import logger
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import jensenshannon, squareform

from prostate_cancer.utils import load_exported_cells, resolve_export_dir

N_PATIENT_CLUSTERS = 6


def core_level_proportions(cells: pd.DataFrame) -> pd.DataFrame:
    """Cell-type proportions per tma_id, pseudocount 1 before normalizing."""
    counts = cells.groupby(["tma_id", "label"], observed=True).size().unstack(fill_value=0)
    counts = counts + 1  # pseudocount, per the paper's Methods
    return counts.div(counts.sum(axis=1), axis=0)


def patient_level_composition(core_proportions: pd.DataFrame, core_to_patient: pd.Series) -> pd.DataFrame:
    """Max-pool core-level proportions to one composition vector per patient."""
    return core_proportions.groupby(core_to_patient).max()


def jsd_linkage(composition: pd.DataFrame) -> np.ndarray:
    n = len(composition)
    dist = np.zeros((n, n))
    values = composition.values
    for i in range(n):
        for j in range(i + 1, n):
            dist[i, j] = dist[j, i] = jensenshannon(values[i], values[j], base=2)
    return linkage(squareform(dist, checks=False), method="average")


def main(export_dir: Path | None = None):
    export_dir = export_dir or resolve_export_dir()
    save_dir = export_dir / "figures" / "figure4"
    save_dir.mkdir(parents=True, exist_ok=True)

    logger.info("loading exported tables")
    cells = load_exported_cells(export_dir, exclude_undefined=True)
    clinical = pd.read_parquet(export_dir / "clinical.parquet")
    assert {"tma_id", "pat_id"} <= set(clinical.columns)

    # sample_id (acquisition) -> tma_id (physical core) -> pat_id: see
    # REPRODUCIBILITY.md for why this de-duplication step is required.
    roi_to_core = clinical["tma_id"]
    roi_to_patient = clinical["pat_id"]
    cells = cells.assign(tma_id=cells["sample_id"].map(roi_to_core))
    assert cells["tma_id"].notna().all(), "some cells' sample_id is missing a clinical/tma_id mapping"

    core_proportions = core_level_proportions(cells)
    assert core_proportions.shape[1] == 34, f"expected 34 cell types, got {core_proportions.shape[1]}"

    core_to_patient = clinical.drop_duplicates("tma_id").set_index("tma_id")["pat_id"]
    core_to_patient = core_to_patient.loc[core_proportions.index]
    composition = patient_level_composition(core_proportions, core_to_patient)
    logger.info(f"patient-level composition matrix: {composition.shape[0]} patients x {composition.shape[1]} cell types")

    Z = jsd_linkage(composition)
    clusters = fcluster(Z, t=N_PATIENT_CLUSTERS, criterion="maxclust")
    patient_clusters = pd.Series(
        [f"P{c}" for c in clusters], index=composition.index, name="patient_cluster"
    )
    logger.info(f"cluster sizes:\n{patient_clusters.value_counts().sort_index()}")

    # reset_index() so `pat_id` is an unambiguous plain column for the R
    # survival script to read, rather than relying on pandas' index metadata
    # round-tripping the same way through arrow's R reader.
    composition.reset_index().to_parquet(save_dir / "figure4a_patient_composition.parquet")
    patient_clusters.reset_index().to_parquet(save_dir / "figure4a_patient_clusters.parquet")

    # %% Figure 4a: clustered heatmap, rows ordered/colored by the 6 patient clusters
    cluster_colors = dict(zip(sorted(patient_clusters.unique()), sns.color_palette("tab10", N_PATIENT_CLUSTERS)))
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
    cg.ax_heatmap.set_title("Figure 4a -- patient-level cell-type composition (P1-P6)")
    cg.figure.savefig(save_dir / "figure4a_composition_heatmap.png", dpi=200)

    logger.info(f"saved figure 4a panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
