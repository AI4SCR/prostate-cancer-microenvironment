"""
This script fits a fresh UMAP, which is NOT reproducible against the
published Figure 2b -- `UMAP.fit()` is never seeded anywhere in the legacy
pipeline, so each fit gives a geometrically different (if topologically
similar) embedding. The actual ground-truth embedding for Figure 2b was
ported once from the legacy `reducer.pkl` and is staged at
`DATA_DIR/umap/all_cells.parquet` -- copy it to
`OUTPUT_FIGURES_DIR/figure2/reducer_embedding.parquet` before running
`figure2_umap.py`, or this script will fit (and cache) a fresh,
non-reproducible embedding instead. See `data/assets.md` and
REPRODUCIBILITY.md for how the ported embedding was produced.

(params match this script's own N_NEIGHBORS/MIN_DIST/NON_MARKER_CHANNELS,
confirmed against `archive/scripts/02-umaps/0-umaps.py`'s `params` list.)
"""

from pathlib import Path

import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import NON_MARKER_CHANNELS, resolve_data_dir, resolve_output_figures_dir

N_NEIGHBORS = 50
MIN_DIST = 0.1


def main(data_dir: Path | None = None):
    import umap

    data_dir = data_dir or resolve_data_dir()
    save_dir = resolve_output_figures_dir() / "figure2"
    save_dir.mkdir(parents=True, exist_ok=True)
    reducer_path = save_dir / "reducer_embedding.parquet"

    if reducer_path.exists():
        logger.info(f"embedding already cached at {reducer_path}, nothing to do")
        return

    logger.info("loading exported tables")
    data = pd.read_parquet(data_dir / "cells" / "intensity_normalized.parquet")
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
