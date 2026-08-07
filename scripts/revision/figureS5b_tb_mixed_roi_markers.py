"""Reviewer response (Issue 5): does nuclear segmentation let signal from
neighboring cells leak into a cell's marker profile? Selects the 3 TMA-core
ROIs with the highest proportion of the ambiguous `immune-T-helper-B-cells`
cell-type label (a population that, if segmentation let T- and B-cell signal
cross-contaminate, would show up exactly as cells scoring positive for both
CD3 and CD20), and renders each ROI as a DNA/CD3/CD20 RGB composite with
nuclear segmentation boundaries overlaid, colored by each cell's own
cell-type call (T cell / B cell / T-B mixed / other), so any signal
spillover across segmentation boundaries is visible directly. DNA is left
out of the composite (for now) since it fills nearly the whole frame and
buries the segmentation outlines against it.

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
from skimage.morphology import dilation, disk
from skimage.segmentation import find_boundaries

from prostate_cancer.utils import normalize_img, resolve_base_dir, resolve_export_dir, resolve_output_figures_dir

TB_MIXED_LABEL = "immune-T-helper-B-cells"
T_LABELS = {
    "immune-T-cells(CD3+)",
    "immune-T-cells_cytotoxic(CD3+CD8a+)",
    "immune-T-cells_helper(CD3+CD4+)",
    "immune-T-cells_regulatory(CD3+CD4+FoxP3+)",
}
B_LABELS = {"immune-B-cells(CD20+)"}
N_ROIS = 3
# R/G channel assignment; no blue/DNA channel for now (see module docstring)
RGB_TARGETS = {"red": "cd20", "green": "cd3"}
OUTLINE_WIDTH = 1  # pixels, via binary dilation of the 1px find_boundaries mask

# Mask-outline color per cell-type category (category code -> RGB, alpha
# fixed at draw time). OTHER covers every label not in T_LABELS/B_LABELS/
# TB_MIXED_LABEL -- light gray at moderate alpha so the full segmentation
# grid stays visible. T/B use soft pastel hues, deliberately distinct from
# the saturated red/green marker channels; TB_MIXED is the population this
# whole figure is about, so it stays bright yellow to stand out rather than
# match the pastel scheme.
OTHER, T, B, TB_MIXED = 0, 1, 2, 3
CATEGORY_COLORS = {
    OTHER: (0.85, 0.85, 0.85),
    T: (0.55, 0.65, 0.95),  # pastel blue
    B: (0.95, 0.55, 0.75),  # pastel pink
    TB_MIXED: (1.0, 1.0, 0.0),  # bright yellow
}
CATEGORY_ALPHA = {OTHER: 0.55, T: 1.0, B: 1.0, TB_MIXED: 1.0}
CATEGORY_LABELS = {OTHER: "other", T: "T cell", B: "B cell", TB_MIXED: "T/B mixed"}

# RGB channel -> marker name, for the marker-color legend
CHANNEL_LABELS = {"red": "CD20", "green": "CD3"}


def categorize_label(label: str) -> int:
    if label == TB_MIXED_LABEL:
        return TB_MIXED
    if label in T_LABELS:
        return T
    if label in B_LABELS:
        return B
    return OTHER


def load_cells(export_dir: Path):
    # cell_annotation.parquet's embedded pandas column-index metadata is
    # stale (predates the sample_id/object_id columns), which makes plain
    # `pd.read_parquet` silently drop them -- ignore_metadata=True forces a
    # column-name-based reconstruction instead.
    cells = pq.read_table(export_dir / "cell_annotation.parquet").to_pandas(ignore_metadata=True)
    return cells[["sample_id", "object_id", "label"]]


def rank_rois_by_tb_mixed_proportion(cells) -> list[str]:
    n_cells = cells.groupby("sample_id").size()
    n_tb_mixed = cells[cells["label"] == TB_MIXED_LABEL].groupby("sample_id").size()
    proportion = (n_tb_mixed / n_cells.reindex(n_tb_mixed.index)).sort_values(ascending=False)

    top_rois = proportion.head(N_ROIS)
    logger.info(f"Top {N_ROIS} ROIs by {TB_MIXED_LABEL} proportion:\n{top_rois}")
    return top_rois.index.tolist()


def object_category_map(cells, sample_id: str, max_object_id: int) -> np.ndarray:
    # lookup[object_id] = category code, indexable directly by mask values;
    # index 0 (background, no cell) stays OTHER/gray but is never a
    # boundary pixel so it's never drawn
    sample_cells = cells[cells["sample_id"] == sample_id]
    lookup = np.full(max_object_id + 1, OTHER, dtype=np.int8)
    categories = sample_cells["label"].map(categorize_label).to_numpy()
    lookup[sample_cells["object_id"].to_numpy()] = categories
    return lookup


def compose_rgb(image: np.ndarray, panel, targets: dict[str, str]) -> np.ndarray:
    # normalize_img arcsinh-transforms and clips at the 99.9th percentile,
    # but doesn't rescale to [0, 1] -- the clip threshold is data-dependent,
    # so without rescaling a channel can come out uniformly dim (threshold
    # < 1) or lose contrast among its brightest cells (threshold > 1, all
    # flattened to 1 by the final clip). Divide by the post-clip max (which
    # equals the clip threshold) so every channel spans the full range.
    target_to_page = panel.reset_index().set_index("target")["channel_index"]
    shape = image.shape[1:]
    channels = {"red": np.zeros(shape), "green": np.zeros(shape), "blue": np.zeros(shape)}
    for color, target in targets.items():
        page = target_to_page.loc[target]
        channel_img = normalize_img(image[page][None, ...])[0]
        channel_max = channel_img.max()
        channels[color] = channel_img / channel_max if channel_max > 0 else channel_img
    rgb = np.stack([channels["red"], channels["green"], channels["blue"]], axis=-1)
    return np.clip(rgb, 0, 1)


def plot_roi(sample_id: str, rgb: np.ndarray, mask: np.ndarray, category_lookup: np.ndarray, save_path: Path):
    boundaries = find_boundaries(mask, mode="inner")
    pixel_category = category_lookup[mask]

    overlay = np.zeros((*boundaries.shape, 4))
    footprint = disk(OUTLINE_WIDTH // 2) if OUTLINE_WIDTH > 1 else None
    # draw OTHER first so thicker T/B/T-B outlines (drawn after) aren't
    # partly overwritten by a neighboring OTHER cell's dilated boundary
    for category in sorted(CATEGORY_COLORS, key=lambda c: c == OTHER, reverse=True):
        color = CATEGORY_COLORS[category]
        pixels = boundaries & (pixel_category == category)
        if footprint is not None:
            pixels = dilation(pixels, footprint)
        overlay[pixels] = (*color, CATEGORY_ALPHA[category])

    # size the figure to the ROI's own aspect ratio (not a fixed square) and
    # render at a high enough DPI that individual nuclei/outlines stay sharp
    h, w = rgb.shape[:2]
    long_edge_in = 10
    figsize = (long_edge_in, long_edge_in * h / w) if w >= h else (long_edge_in * w / h, long_edge_in)

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(rgb, interpolation="none")
    ax.imshow(overlay, interpolation="none")

    marker_handles = [
        plt.Line2D([0], [0], color=color, lw=4, label=CHANNEL_LABELS[color_name])
        for color_name, color in {"red": (1, 0, 0), "green": (0, 1, 0)}.items()
    ]
    marker_legend = ax.legend(handles=marker_handles, loc="upper left", fontsize=8, framealpha=0.6, title="Marker", title_fontsize=8)
    ax.add_artist(marker_legend)

    category_handles = [
        plt.Line2D([0], [0], color=color, lw=3, label=CATEGORY_LABELS[category]) for category, color in CATEGORY_COLORS.items()
    ]
    ax.legend(handles=category_handles, loc="upper right", fontsize=8, framealpha=0.6, title="Nuclear outline", title_fontsize=8)

    ax.set_title(sample_id)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(save_path, dpi=400)
    plt.close(fig)


def main(base_dir: Path | None = None, export_dir: Path | None = None):
    base_dir = base_dir or resolve_base_dir()
    export_dir = export_dir or resolve_export_dir()
    save_dir = resolve_output_figures_dir().parent / "revision" / "figureS5b_tb_mixed_roi_markers"
    save_dir.mkdir(parents=True, exist_ok=True)

    cells = load_cells(export_dir)
    sample_ids = rank_rois_by_tb_mixed_proportion(cells)

    dataset = PCa(base_dir=base_dir, image_version="filtered", mask_version="annotated", load_intensity=False, load_metadata=False, align=False)
    dataset.setup(sample_ids=sample_ids)

    for sample_id in sample_ids:
        image = dataset.images[sample_id].data
        mask = dataset.masks[sample_id].data
        rgb = compose_rgb(image, dataset.panel, RGB_TARGETS)
        category_lookup = object_category_map(cells, sample_id, max_object_id=int(mask.max()))
        plot_roi(sample_id, rgb, mask, category_lookup, save_dir / f"{sample_id}.png")

    logger.info(f"Saved {len(sample_ids)} ROI marker panels to {save_dir}")


if __name__ == "__main__":
    CLI(main)
