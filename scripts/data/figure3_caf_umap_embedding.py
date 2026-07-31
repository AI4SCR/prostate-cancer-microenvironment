# %%
"""Compute and cache the two UMAP embeddings `scripts/figures/figure3_caf_umap.py` plots.

Split out of the original combined fit+plot script so the (fast) plotting
script can be re-run for styling changes without re-triggering these slow
UMAP fits. Pure data-loading/export mechanics reorganization -- see
`scripts/figures/figure3_caf_umap.py`'s docstring for the full port
provenance (cell filter, per-config normalization, both marker configs); the
fits themselves are unchanged.

Reads from $EXPORT_DIR. Writes each config's embedding to
$EXPORT_DIR/figures/figure3/{config_name}/umap_embedding.parquet.
"""
from pathlib import Path

import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import NON_MARKER_CHANNELS, normalize, resolve_export_dir

ALL_MARKERS = [
    "smooth_muscle_actin", "prostate_specific_antigen", "vimentin", "collagen1", "synaptophysin", "keratin5",
    "yap1", "pan_keratin", "ces1", "egr1", "cd31", "cd45", "cd44", "fap", "fox_p3", "cd4", "e_cadherin", "cd68",
    "cd66b", "cd20", "cd8a", "cd11b", "p63", "beta_catenin", "pdpn", "cd105", "ki_67", "p53", "cd3", "erg",
    "c_casp3", "cnn1", "keratin8_18", "cd146", "ar", "dna1", "dna2", "icsk1", "icsk2", "icsk3",
]
CAF_MARKERS = {"smooth_muscle_actin", "vimentin", "collagen1", "cd146", "cnn1", "cd105", "ar", "egr1", "ces1"}
CONFIGS = {
    "excl_markers": sorted(NON_MARKER_CHANNELS),
    "caf_markers_only": sorted(set(ALL_MARKERS) - CAF_MARKERS),
}
N_NEIGHBORS = 50
MIN_DIST = 0.1


def main(export_dir: Path | None = None):
    import umap

    export_dir = export_dir or resolve_export_dir()
    save_dir = export_dir / "figures" / "figure3"
    save_dir.mkdir(parents=True, exist_ok=True)

    metadata = None
    raw_intensity = None

    for config_name, exclude_markers in CONFIGS.items():
        config_dir = save_dir / config_name
        config_dir.mkdir(parents=True, exist_ok=True)
        embedding_path = config_dir / "umap_embedding.parquet"

        if embedding_path.exists():
            logger.info(f"[{config_name}] embedding already cached at {embedding_path}, nothing to do")
            continue

        if metadata is None:
            logger.info("loading exported tables")
            metadata = pd.read_parquet(export_dir / "metadata.parquet")
            raw_intensity = pd.read_parquet(export_dir / "intensity.parquet")
            metadata, raw_intensity = metadata.align(raw_intensity, axis=0, join="inner")
            caf_filter = metadata["label"].str.contains("CAF")
            raw_intensity = raw_intensity.loc[caf_filter, :]
            logger.info(f"{len(raw_intensity)} CAF cells (label contains 'CAF')")

        df = raw_intensity.loc[:, ~raw_intensity.columns.isin(exclude_markers)].copy()
        df = normalize(df, exclude_zeros=True)

        logger.info(f"[{config_name}] computing UMAP on {len(df)} CAF cells x {df.shape[1]} markers")
        reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, metric="euclidean")
        reducer.fit(df.values)
        embedding_df = pd.DataFrame(reducer.embedding_, index=df.index, columns=["umap_1", "umap_2"])
        embedding_df.to_parquet(embedding_path)
        logger.info(f"[{config_name}] saved embedding to {embedding_path}")


if __name__ == "__main__":
    CLI(main)
