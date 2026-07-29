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

import matplotlib.pyplot as plt
import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import create_color_maps, resolve_export_dir

# Non-biological channels present in intensity_normalized.parquet: DNA
# intercalator + cell-segmentation-kit channels, plus FAP (excluded by the
# paper's Methods due to non-specific staining after in-house conjugation).
NON_MARKER_COLUMNS = ["dna1", "dna2", "icsk1", "icsk2", "icsk3", "fap"]
INDEX_COLUMNS = ["sample_id", "object_id", "slide_code", "donor_block_id", "pat_id"]
LABEL_COLUMNS = ["label", "main_group", "label_id", "main_group_id", "meta_label", "meta_label_id"]


def load_cells(export_dir: Path) -> pd.DataFrame:
    metadata = pd.read_parquet(export_dir / "metadata.parquet").reset_index()
    intensity = pd.read_parquet(export_dir / "intensity_normalized.parquet").reset_index()
    cells = intensity.merge(metadata, on=["sample_id", "object_id"], validate="one_to_one")
    return cells


def plot_categorical(cells: pd.DataFrame, col: str, save_path: Path, point_size: float = 2.0):
    color_maps = create_color_maps(cells[[col]].astype({col: "category"}))
    color_map = color_maps[col]

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(
        cells["umap_1"], cells["umap_2"],
        c=[color_map[v] for v in cells[col]],
        s=point_size, alpha=0.5, linewidths=0,
    )
    ax.set_title(f"Figure 2b -- UMAP colored by {col}")
    ax.set_xticks([])
    ax.set_yticks([])

    # too many categories (e.g. 195 patients) to legend legibly -- skip it
    if len(color_map) <= 40:
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="", color=c, label=str(l))
            for l, c in color_map.items()
        ]
        ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6, ncols=2 if len(handles) > 15 else 1)

    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def plot_marker(cells: pd.DataFrame, marker: str, save_path: Path, point_size: float = 2.0):
    fig, ax = plt.subplots(figsize=(8, 7))
    sca = ax.scatter(cells["umap_1"], cells["umap_2"], c=cells[marker], s=point_size, alpha=0.6, cmap="viridis", linewidths=0)
    fig.colorbar(sca, ax=ax, label=marker)
    ax.set_title(f"Figure 2b -- UMAP colored by {marker}")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


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
    cells = load_cells(export_dir)

    marker_cols = [c for c in cells.columns if c not in NON_MARKER_COLUMNS + INDEX_COLUMNS + LABEL_COLUMNS]
    assert len(marker_cols) == 34, f"expected 34 markers, got {len(marker_cols)}: {marker_cols}"
    assert set(markers) <= set(marker_cols), f"{set(markers) - set(marker_cols)} not in marker panel"

    if 0 < n_cells < len(cells):
        logger.info(f"subsampling {n_cells} of {len(cells)} cells (seed={seed})")
        cells = cells.sample(n=n_cells, random_state=seed).reset_index(drop=True)

    logger.info(f"computing UMAP on {len(cells)} cells x {len(marker_cols)} markers")
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric="euclidean", random_state=seed)
    embedding = reducer.fit_transform(cells[marker_cols].values)
    cells["umap_1"] = embedding[:, 0]
    cells["umap_2"] = embedding[:, 1]

    cells.to_parquet(save_dir / "umap_embedding.parquet")

    for col in ["main_group", "label", "pat_id"]:
        plot_categorical(cells, col, save_dir / f"umap_{col}.png")

    for marker in markers:
        plot_marker(cells, marker, save_dir / f"umap_marker_{marker}.png")

    logger.info(f"saved figure 2b panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
