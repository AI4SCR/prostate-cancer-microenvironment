"""
This script fits a fresh UMAP, which is NOT reproducible against the
published Figure 3a -- `UMAP.fit()` is never seeded anywhere in the legacy
pipeline, so each fit gives a geometrically different (if topologically
similar) embedding. The actual ground-truth embeddings for Figure 3a were
ported once from the legacy `reducer.pkl` files and are staged at
`DATA_DIR/umap/caf_{excl_markers,markers_only,stromal}.parquet` -- copy the
relevant one to `OUTPUT_FIGURES_DIR/figure3/<config>/umap_embedding.parquet`
before running `figure3_caf_umap.py`, or this script will fit (and cache) a
fresh, non-reproducible embedding instead. See `data/assets.md` and
REPRODUCIBILITY.md for how the ported embeddings were produced.

(params match this script's own N_NEIGHBORS/MIN_DIST/CONFIGS, confirmed
against `archive/scripts/02-umaps/0-umaps-cafs.py`'s `params` list -- the
legacy dir names truncate the marker list to 150 chars, hence the shorter
`caf_markers_only` name above despite excluding more markers.

The published Figure 3a legend text says "UMAP of all CAF cells" -- this is
a manuscript error, confirmed directly: the panel shown is actually the
full stromal-compartment UMAP (`2-main_groups/main_group=stromal`), which
includes `stromal-pericytes` alongside the CAF subtypes. The "stromal"
config below ports that reducer; `excl_markers`/`caf_markers_only` above
are the CAF-only reducers, kept for reference/other panels but no longer
the source of the main Figure 3a UMAP.)
"""

from pathlib import Path

import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import NON_MARKER_CHANNELS, normalize, resolve_data_dir, resolve_output_figures_dir

ALL_MARKERS = [
    "smooth_muscle_actin", "prostate_specific_antigen", "vimentin", "collagen1", "synaptophysin", "keratin5",
    "yap1", "pan_keratin", "ces1", "egr1", "cd31", "cd45", "cd44", "fap", "fox_p3", "cd4", "e_cadherin", "cd68",
    "cd66b", "cd20", "cd8a", "cd11b", "p63", "beta_catenin", "pdpn", "cd105", "ki_67", "p53", "cd3", "erg",
    "c_casp3", "cnn1", "keratin8_18", "cd146", "ar", "dna1", "dna2", "icsk1", "icsk2", "icsk3",
]
CAF_MARKERS = {"smooth_muscle_actin", "vimentin", "collagen1", "cd146", "cnn1", "cd105", "ar", "egr1", "ces1"}
# (exclude_markers, population) per config -- population is "caf" (label contains "CAF")
# or "stromal" (main_group == "stromal", includes pericytes and other non-CAF stromal cells)
CONFIGS = {
    "excl_markers": (sorted(NON_MARKER_CHANNELS), "caf"),
    "caf_markers_only": (sorted(set(ALL_MARKERS) - CAF_MARKERS), "caf"),
    "stromal": (sorted(NON_MARKER_CHANNELS), "stromal"),
}
N_NEIGHBORS = 50
MIN_DIST = 0.1


def main(data_dir: Path | None = None):
    import umap

    data_dir = data_dir or resolve_data_dir()
    save_dir = resolve_output_figures_dir() / "figure3"
    save_dir.mkdir(parents=True, exist_ok=True)

    metadata = None
    raw_intensity = None

    for config_name, (exclude_markers, population) in CONFIGS.items():
        config_dir = save_dir / config_name
        config_dir.mkdir(parents=True, exist_ok=True)
        embedding_path = config_dir / "umap_embedding.parquet"

        if embedding_path.exists():
            logger.info(f"[{config_name}] embedding already cached at {embedding_path}, nothing to do")
            continue

        if metadata is None:
            logger.info("loading exported tables")
            metadata = pd.read_parquet(data_dir / "cells" / "metadata.parquet")
            raw_intensity = pd.read_parquet(data_dir / "cells" / "intensity.parquet")
            metadata, raw_intensity = metadata.align(raw_intensity, axis=0, join="inner")

        population_filter = metadata["label"].str.contains("CAF") if population == "caf" else metadata["main_group"] == "stromal"
        raw_intensity_subset = raw_intensity.loc[population_filter, :]
        logger.info(f"[{config_name}] {len(raw_intensity_subset)} {population} cells")

        df = raw_intensity_subset.loc[:, ~raw_intensity_subset.columns.isin(exclude_markers)].copy()
        df = normalize(df, exclude_zeros=True)

        logger.info(f"[{config_name}] computing UMAP on {len(df)} {population} cells x {df.shape[1]} markers")
        reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, metric="euclidean")
        reducer.fit(df.values)
        embedding_df = pd.DataFrame(reducer.embedding_, index=df.index, columns=["umap_1", "umap_2"])
        embedding_df.to_parquet(embedding_path)
        logger.info(f"[{config_name}] saved embedding to {embedding_path}")


if __name__ == "__main__":
    CLI(main)
