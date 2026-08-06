"""
This script fits a fresh UMAP, which is NOT reproducible against the
published Figure 2b -- `UMAP.fit()` is never seeded anywhere in the legacy
pipeline, so each fit gives a geometrically different (if topologically
similar) embedding. The actual ground-truth embedding for Figure 2b is
ported from the legacy `reducer.pkl` via `scripts/port/port_umap_reducer.py`
(run from that script's own pinned pixi env, NOT this repo's `.venv`):

    cd scripts/port
    pixi run python port_umap_reducer.py \\
      "/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/2-umaps/0-all-cells/n_neighbors=50-min_dist=0.1-engine=umap-learn-excl_markers=dna1_dna2_icsk1_icsk2_icsk3_fap/reducer.pkl" \\
      ../../data/figures/figure2_umap/umap_embeddings.parquet

(params match this script's own N_NEIGHBORS/MIN_DIST/NON_MARKER_CHANNELS,
confirmed against `archive/scripts/02-umaps/0-umaps.py`'s `params` list.)
"""

from pathlib import Path

import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import NON_MARKER_CHANNELS, resolve_export_dir, resolve_output_figures_dir

N_NEIGHBORS = 50
MIN_DIST = 0.1


def main(export_dir: Path | None = None):
    import umap

    export_dir = export_dir or resolve_export_dir()
    save_dir = resolve_output_figures_dir() / "figure2"
    save_dir.mkdir(parents=True, exist_ok=True)
    reducer_path = save_dir / "reducer_embedding.parquet"

    if reducer_path.exists():
        logger.info(f"embedding already cached at {reducer_path}, nothing to do")
        return

    logger.info("loading exported tables")
    data = pd.read_parquet(export_dir / "intensity_normalized.parquet")
    fit_data = data.loc[:, ~data.columns.isin(NON_MARKER_CHANNELS)]

    logger.info(
        f"computing UMAP for {len(fit_data)} cells, n_neighbors={N_NEIGHBORS}, "
        f"min_dist={MIN_DIST}, excluding {NON_MARKER_CHANNELS}"
    )
    reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, metric="euclidean")
    reducer.fit(fit_data.values)
    embedding_df = pd.DataFrame(reducer.embedding_, index=fit_data.index, columns=["umap_1", "umap_2"])
    embedding_df.to_parquet(reducer_path)
    logger.info(f"saved embedding to {reducer_path}")


if __name__ == "__main__":
    CLI(main)
