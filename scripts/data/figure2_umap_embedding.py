# %%
"""Compute and cache the UMAP embedding `scripts/figures/figure2_umap.py` plots.

Split out of the original combined fit+plot script so the (fast) plotting
script can be re-run for styling changes without re-triggering this slow
(~2.19M-cell) `UMAP.fit()`. Pure data-loading/export mechanics reorganization
-- see `scripts/figures/figure2_umap.py`'s docstring for the full port
provenance; the fit itself (parameters, data, seeding behavior) is unchanged.

Reads from $EXPORT_DIR. Writes the fitted embedding to
$EXPORT_DIR/figures/figure2/reducer_embedding.parquet.
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
    # No random_state -- legacy's compute_umap()/run_umap() never passes one
    # to UMAP() for the fit itself (only figure2_umap.py's later
    # subsampling-for-plotting step is seeded). Also lets UMAP run
    # multi-threaded instead of the single-threaded path a fixed
    # random_state forces.
    reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, metric="euclidean")
    reducer.fit(fit_data.values)
    embedding_df = pd.DataFrame(reducer.embedding_, index=fit_data.index, columns=["umap_1", "umap_2"])
    embedding_df.to_parquet(reducer_path)
    logger.info(f"saved embedding to {reducer_path}")


if __name__ == "__main__":
    CLI(main)
