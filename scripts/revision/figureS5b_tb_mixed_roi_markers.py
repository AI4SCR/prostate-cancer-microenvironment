"""Reviewer response (Issue 5): does nuclear segmentation let signal from
neighboring cells leak into a cell's marker profile? Selects the 3 TMA-core
ROIs with the highest proportion of the ambiguous `immune-T-helper-B-cells`
cell-type label (a population that, if segmentation let T- and B-cell signal
cross-contaminate, would show up exactly as cells scoring positive for both
CD3 and CD20), and renders each ROI as a DNA/CD3/CD20 RGB composite with the
nuclear segmentation boundaries overlaid, so any signal spillover across
segmentation boundaries is visible directly.

New figure (no legacy port): the cell-type label used to rank ROIs already
exists in `cell_annotation.parquet` (produced by `export.py`), and the raw
images/masks read here are staged, unmodified dataset files
(`02_processed/images/filtered`, `02_processed/masks/annotated` under
`BASE_DIR`) -- nothing here derives new biological data, only selects and
renders what already exists.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
from ai4bmr_datasets import PCa
from jsonargparse import CLI
from loguru import logger
from skimage.segmentation import find_boundaries

from prostate_cancer.utils import normalize_img, resolve_base_dir, resolve_export_dir, resolve_output_figures_dir

TB_MIXED_LABEL = "immune-T-helper-B-cells"
N_ROIS = 3
# R/G/B channel assignment: standard IF convention (blue nuclei, marker
# channels in red/green)
RGB_TARGETS = {"red": "cd20", "green": "cd3", "blue": "dna1"}


def rank_rois_by_tb_mixed_proportion(export_dir: Path) -> list[str]:
    # cell_annotation.parquet's embedded pandas column-index metadata is
    # stale (predates the sample_id/object_id columns), which makes plain
    # `pd.read_parquet` silently drop them -- ignore_metadata=True forces a
    # column-name-based reconstruction instead.
    cells = pq.read_table(export_dir / "cell_annotation.parquet").to_pandas(ignore_metadata=True)
    cells = cells[["sample_id", "label"]]

    n_cells = cells.groupby("sample_id").size()
    n_tb_mixed = cells[cells["label"] == TB_MIXED_LABEL].groupby("sample_id").size()
    proportion = (n_tb_mixed / n_cells.reindex(n_tb_mixed.index)).sort_values(ascending=False)

    top_rois = proportion.head(N_ROIS)
    logger.info(f"Top {N_ROIS} ROIs by {TB_MIXED_LABEL} proportion:\n{top_rois}")
    return top_rois.index.tolist()


def compose_rgb(image: np.ndarray, panel, targets: dict[str, str]) -> np.ndarray:
    target_to_page = panel.reset_index().set_index("target")["channel_index"]
    channels = {}
    for color, target in targets.items():
        page = target_to_page.loc[target]
        channel_img = normalize_img(image[page][None, ...])[0]
        channels[color] = channel_img
    rgb = np.stack([channels["red"], channels["green"], channels["blue"]], axis=-1)
    return np.clip(rgb, 0, 1)


def plot_roi(sample_id: str, rgb: np.ndarray, mask: np.ndarray, save_path: Path):
    boundaries = find_boundaries(mask, mode="inner")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(rgb)
    overlay = np.zeros((*boundaries.shape, 4))
    overlay[boundaries] = [1, 1, 1, 0.6]
    ax.imshow(overlay)
    ax.set_title(f"{sample_id}\nR=CD20 G=CD3 B=DNA1, nuclear boundaries overlaid")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)


def main(base_dir: Path | None = None, export_dir: Path | None = None):
    base_dir = base_dir or resolve_base_dir()
    export_dir = export_dir or resolve_export_dir()
    save_dir = resolve_output_figures_dir().parent / "revision" / "figureS5b_tb_mixed_roi_markers"
    save_dir.mkdir(parents=True, exist_ok=True)

    sample_ids = rank_rois_by_tb_mixed_proportion(export_dir)

    dataset = PCa(base_dir=base_dir, image_version="filtered", mask_version="annotated", load_intensity=False, load_metadata=False, align=False)
    dataset.setup(sample_ids=sample_ids)

    for sample_id in sample_ids:
        image = dataset.images[sample_id].data
        mask = dataset.masks[sample_id].data
        rgb = compose_rgb(image, dataset.panel, RGB_TARGETS)
        plot_roi(sample_id, rgb, mask, save_dir / f"{sample_id}.png")

    logger.info(f"Saved {len(sample_ids)} ROI marker panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
