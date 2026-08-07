"""Reproduce Supplementary Figure 7: 50-seed k-means ARI robustness sweep.

1:1 port of `sync_paper/06-spatial-niches/construction/kmeans_clustering.py`'s
ARI-robustness section (lines ~217-282 of that file) -- the same file whose
final section (`k=24, seed=686`) is already ported as
`figure5_niche_clustering.py`'s single clustering step. This script instead
runs k-means once per each of 50 seeds (`np.random.RandomState(42)`),
computes the pairwise Adjusted Rand Index across all 50 runs, and boxplots
ARI per run -- the best-agreement run highlighted yellow, which is
literally `seed=686`, the same seed `figure5_niche_clustering.py` uses for
the actual published clustering (confirmed in the legacy source's own
comment: "seed with top ARI mean across runs").

Does not touch `figure5_niche_clustering.py` -- separate, dedicated script,
same data/utils, per this project's established pattern for panels sharing
a legacy file with a different figure's clustering step.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from jsonargparse import CLI
from loguru import logger
from matplotlib import pyplot as plt

from prostate_cancer.utils import resolve_data_dir, resolve_output_figures_dir

K = 24
N_RUNS = 50
SEED_RNG = 42
GRAPH_RADIUS = 32
NEIGHBOR_COUNT_THRESHOLD = 2  # drop cells with <= this many neighbors
MIN_CELLS_PER_NICHE = 15
MIN_PATIENTS_PER_NICHE = 5


def construct_ari_matrix(df_results):
    ari_matrix = pd.DataFrame(index=df_results.columns, columns=df_results.columns, dtype=float)
    from sklearn.metrics import adjusted_rand_score

    for i, col1 in enumerate(df_results.columns):
        for j, col2 in enumerate(df_results.columns):
            if i == j:
                ari_matrix.loc[col1, col2] = 1.0
            elif pd.isna(ari_matrix.loc[col1, col2]):
                ari = adjusted_rand_score(df_results[col1], df_results[col2])
                ari_matrix.loc[col1, col2] = ari
                ari_matrix.loc[col2, col1] = ari
    return ari_matrix


def prepare_ari_data(ari_matrix):
    ari_results = pd.DataFrame(index=ari_matrix.index, columns=[f"run{i+1}" for i in range(len(ari_matrix) - 1)], dtype=float)
    for i, col in enumerate(ari_matrix.columns):
        entries = ari_matrix[col].values
        entries = np.delete(entries, i)
        ari_results.iloc[i] = entries
    ari_results = ari_results.T
    ari_melted = ari_results.melt(var_name="clustering", value_name="ARI")
    return ari_melted


def main(data_dir: Path | None = None):
    data_dir = data_dir or resolve_data_dir()
    neighborhoods_dir = data_dir / "neighborhoods"
    save_dir = resolve_output_figures_dir() / "figureS7"
    save_dir.mkdir(parents=True, exist_ok=True)

    sys.path.append(str(data_dir / "niches" / "robustness"))
    from utils.clustering import perform_kmeans_clustering, wrapper_nhood_filtering

    # %% load the neighborhood-composition graph (same data as figure5_niche_clustering.py)
    data_path = neighborhoods_dir / f"radius{GRAPH_RADIUS}_data.parquet"
    cell_metadata = pd.read_parquet(neighborhoods_dir / "cell_metadata.parquet", engine="fastparquet")
    metadata = pd.read_parquet(neighborhoods_dir / "metadata.parquet", engine="fastparquet")

    data = pd.read_parquet(data_path, engine="fastparquet")
    data_filt = data[data.sum(axis=1) > NEIGHBOR_COUNT_THRESHOLD]
    logger.info(f"filtered out {data.shape[0] - data_filt.shape[0]} cells with <= {NEIGHBOR_COUNT_THRESHOLD} neighbors")
    data = data_filt.div(data_filt.sum(axis=1), axis=0)  # counts -> frequencies

    # %% 50-seed k-means sweep
    seeds = np.random.RandomState(SEED_RNG).randint(0, 1000, size=N_RUNS)
    df_results = pd.DataFrame(index=data.index)

    for seed in seeds:
        logger.info(f"running k-means with k={K}, seed={seed}")
        cluster_data, cluster_name = perform_kmeans_clustering(data, k=K, random_state=seed)
        cluster_data = wrapper_nhood_filtering(
            cluster_data, cluster_name, metadata, cell_metadata, min_cells=MIN_CELLS_PER_NICHE, min_patients=MIN_PATIENTS_PER_NICHE
        )
        final_cluster_name = cluster_data.columns[-1]
        df_results = df_results.join(cluster_data[[final_cluster_name]], how="inner")
        assert data.index.equals(df_results.index), f"data index {data.index} does not match df_results index {df_results.index}"
        assert data.shape[0] == df_results.shape[0], f"data shape {data.shape} does not match df_results shape {df_results.shape}"
        logger.info(f"df_results shape: {df_results.shape}")

    # %% robustness checks
    ari_matrix = construct_ari_matrix(df_results)
    ari_results = prepare_ari_data(ari_matrix)
    ari_means = ari_results.groupby("clustering")["ARI"].mean()
    top_run = ari_means.idxmax()

    ari_melted = ari_results.copy()
    palette = {run: ("yellow" if run == top_run else "lightblue") for run in ari_melted["clustering"].unique()}

    plt.figure(figsize=(25, 8))
    ax = plt.gca()
    sns.boxplot(data=ari_melted, x="clustering", y="ARI", palette=palette, ax=ax)

    for i, run in enumerate(ax.get_xticklabels()):
        name = run.get_text()
        v = ari_means[name]
        ax.text(i, v + 0.01, f"{v:.2f}", ha="center", fontsize=10)

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plot_path = save_dir / "ari_boxplot.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    logger.info(f"saved Supplementary Figure 7 to {plot_path}")


if __name__ == "__main__":
    CLI(main)
