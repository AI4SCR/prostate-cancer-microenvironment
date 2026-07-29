# %%
"""Reproduce Figure 2b: UMAP of all cells colored by compartment/cell type/patient/markers.

Reads the tables `scripts/00-data-export/export_for_r.py` produces in
`$EXPORT_DIR` (never `$BASE_DIR` -- see REPRODUCIBILITY.md). Writes PNGs and
the embedding itself to `$EXPORT_DIR/figures/figure2/`.

The paper computed Parametric UMAP via GPU-accelerated rapids_singlecell on
all 2.19M cells. This environment has no GPU, so `n_cells` subsamples for a
CPU-tractable run with umap-learn; pass `n_cells=-1` for the full dataset
(slow without a GPU).
"""
from pathlib import Path

from jsonargparse import CLI
from loguru import logger

from prostate_cancer.plotting import plot_embedding_categorical, plot_embedding_marker
from prostate_cancer.utils import load_exported_cells, marker_columns, resolve_export_dir


def main(
    export_dir: Path | None = None,
    n_cells: int = 50_000,
    seed: int = 0,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    markers: tuple[str, ...] = ("cd45", "pan_keratin", "cd31", "vimentin"),
):
    import umap

    export_dir = export_dir or resolve_export_dir()
    save_dir = export_dir / "figures" / "figure2"
    save_dir.mkdir(parents=True, exist_ok=True)

    logger.info("loading exported tables")
    # Figure 2b's own caption states "n = 2'191'967" -- the full, unfiltered
    # cell count -- so unlike Figure 2a's heatmap, "undefined" is kept.
    cells = load_exported_cells(export_dir, exclude_undefined=False)
    cols = marker_columns(cells)
    assert set(markers) <= set(cols), f"{set(markers) - set(cols)} not in marker panel"

    if 0 < n_cells < len(cells):
        logger.info(f"subsampling {n_cells} of {len(cells)} cells (seed={seed})")
        cells = cells.sample(n=n_cells, random_state=seed).reset_index(drop=True)

    logger.info(f"computing UMAP on {len(cells)} cells x {len(cols)} markers")
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric="euclidean", random_state=seed)
    embedding = reducer.fit_transform(cells[cols].values)
    cells["umap_1"] = embedding[:, 0]
    cells["umap_2"] = embedding[:, 1]

    cells.to_parquet(save_dir / "umap_embedding.parquet")

    for col in ["main_group", "label", "pat_id"]:
        plot_embedding_categorical(cells, col, save_dir / f"umap_{col}.png", title=f"Figure 2b -- UMAP colored by {col}")

    for marker in markers:
        plot_embedding_marker(cells, marker, save_dir / f"umap_marker_{marker}.png", title=f"Figure 2b -- UMAP colored by {marker}")

    logger.info(f"saved figure 2b panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
