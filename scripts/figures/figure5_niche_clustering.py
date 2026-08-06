import sys
from pathlib import Path

import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import resolve_export_dir, resolve_legacy_dir, resolve_output_figures_dir

K = 24
SEED = 686
GRAPH_RADIUS = 32
NEIGHBOR_COUNT_THRESHOLD = 2  # drop cells with <= this many neighbors
MIN_CELLS_PER_NICHE = 15
MIN_PATIENTS_PER_NICHE = 5


def main(export_dir: Path | None = None, legacy_dir: Path | None = None):
    export_dir = export_dir or resolve_export_dir()
    legacy_dir = legacy_dir or resolve_legacy_dir()
    count_base_dir = legacy_dir / "PCa_NHood" / "CellCellNeighborhoods"
    save_dir = resolve_output_figures_dir() / "figure5"
    save_dir.mkdir(parents=True, exist_ok=True)

    sys.path.append(str(legacy_dir / "PCA_NHOODs_clean" / "robustness"))
    from utils.clustering import perform_kmeans_clustering, wrapper_nhood_filtering
    from utils.visualization import calculate_freqs, calculate_zscore, create_plot, transform_for_heatmap

    # %% load the neighborhood-composition graph
    data_path = count_base_dir / f"graph_type=radius-radius={GRAPH_RADIUS}" / "data.parquet"
    cell_metadata = pd.read_parquet(count_base_dir / "cell_metadata.parquet", engine="fastparquet")
    metadata = pd.read_parquet(count_base_dir / "metadata.parquet", engine="fastparquet")

    data = pd.read_parquet(data_path, engine="fastparquet")
    data_filt = data[data.sum(axis=1) > NEIGHBOR_COUNT_THRESHOLD]
    logger.info(f"filtered out {data.shape[0] - data_filt.shape[0]} cells with <= {NEIGHBOR_COUNT_THRESHOLD} neighbors")
    data = data_filt.div(data_filt.sum(axis=1), axis=0)  # counts -> frequencies

    # %% k-means clustering + small-niche filtering
    logger.info(f"running k-means with k={K}, seed={SEED}")
    cluster_data, cluster_name = perform_kmeans_clustering(data, k=K, random_state=SEED)
    cluster_data = wrapper_nhood_filtering(
        cluster_data, cluster_name, metadata, cell_metadata, min_cells=MIN_CELLS_PER_NICHE, min_patients=MIN_PATIENTS_PER_NICHE
    )
    final_cluster_name = cluster_data.columns[-1]

    df_results = pd.DataFrame(index=data.index).join(cluster_data[[final_cluster_name]], how="inner")
    assert data.index.equals(df_results.index), f"data index {data.index} does not match df_results index {df_results.index}"
    assert data.shape[0] == df_results.shape[0], f"data shape {data.shape} does not match df_results shape {df_results.shape}"
    logger.info(f"cluster assignments: {df_results.shape}")

    # %% z-scored niche x cell-type composition heatmap
    df_clusters = df_results.join(cell_metadata, how="right")
    df_clusters.fillna("unassigned", inplace=True)
    df_clusters.reset_index(inplace=True)
    nhood = final_cluster_name
    df_clusters[nhood] = df_clusters[nhood].astype("category")
    df_clusters["label"] = df_clusters["label"].astype("category")

    _, df_freqs = calculate_freqs(df_clusters, cell_type_col="label", nhood_col=nhood, sample_col="sample_id")
    target = "z_score_sample_normalized"
    freq_df = calculate_zscore(
        df=df_freqs, x_col="mean_nhood_celltype_sample", mean_col="mean_celltype_sample", std_col="std_celltype_sample", colname=target
    )
    heatmap_data = transform_for_heatmap(freq_df, target_col=target, cell_type_col="label", nhood_col=nhood)

    p = create_plot(
        df=heatmap_data,
        metric="Frequency normalized sample-wise (z-score)",
        nhood=nhood,
        color_scheme="coolwarm",
        upper_limit=3,
        lower_limit=-3,
        cluster_celltypes=False,
        cluster_neighborhoods=True,
    )
    p.savefig(save_dir / f"{nhood}_heatmap.png", dpi=300, bbox_inches="tight")

    # %%
    df_clusters.set_index(["sample_id", "object_id"], inplace=True)
    df_clusters.to_parquet(save_dir / "clusters.parquet", engine="fastparquet")
    logger.info(f"saved clustering results to {save_dir}")


if __name__ == "__main__":
    CLI(main)
