# %%
"""Reproduce Figure 3a: UMAP of stromal (CAF/SMC) cells colored by subcluster.

Paper: "UMAP of all CAF cells across the cohort, computed using CAF-relevant
markers (alpha-SMA, vimentin, collagen I, CD146, CNN1, CD105, AR, EGR1, and
CES1)." The 10 stromal labels in `metadata.parquet` (all under
main_group=="stromal") are exactly the "10 distinct CAF populations" the
paper's Results describes -- no further filtering needed.

Reads from `$EXPORT_DIR` (never `$BASE_DIR` -- see REPRODUCIBILITY.md).
Writes to `$EXPORT_DIR/figures/figure3/`.
"""
from pathlib import Path

from jsonargparse import CLI
from loguru import logger

from prostate_cancer.plotting import plot_embedding_categorical
from prostate_cancer.utils import load_exported_cells, resolve_export_dir

CAF_MARKERS = ["smooth_muscle_actin", "vimentin", "collagen1", "cd146", "cnn1", "cd105", "ar", "egr1", "ces1"]


def main(
    export_dir: Path | None = None,
    n_cells: int = -1,  # stromal compartment is ~30% of 2.19M cells but still CPU-tractable at full size
    seed: int = 0,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
):
    import umap

    export_dir = export_dir or resolve_export_dir()
    save_dir = export_dir / "figures" / "figure3"
    save_dir.mkdir(parents=True, exist_ok=True)

    logger.info("loading exported tables")
    cells = load_exported_cells(export_dir, exclude_undefined=True)
    cells = cells[cells["main_group"] == "stromal"].reset_index(drop=True)
    n_caf_labels = cells["label"].nunique()
    assert n_caf_labels == 10, f"expected 10 stromal/CAF labels, got {n_caf_labels}"

    if 0 < n_cells < len(cells):
        logger.info(f"subsampling {n_cells} of {len(cells)} stromal cells (seed={seed})")
        cells = cells.sample(n=n_cells, random_state=seed).reset_index(drop=True)

    logger.info(f"computing UMAP on {len(cells)} stromal cells x {len(CAF_MARKERS)} CAF markers")
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric="euclidean", random_state=seed)
    embedding = reducer.fit_transform(cells[CAF_MARKERS].values)
    cells["umap_1"] = embedding[:, 0]
    cells["umap_2"] = embedding[:, 1]

    cells.to_parquet(save_dir / "figure3a_umap_embedding.parquet")
    plot_embedding_categorical(cells, "label", save_dir / "figure3a_umap_label.png", title="Figure 3a -- CAF UMAP colored by subcluster")

    logger.info(f"saved figure 3a panel to {save_dir}")


if __name__ == "__main__":
    CLI(main)
