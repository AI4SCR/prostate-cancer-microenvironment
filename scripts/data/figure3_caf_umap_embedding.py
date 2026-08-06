"""
This script fits a fresh UMAP, which is NOT reproducible against the
published Figure 3a -- `UMAP.fit()` is never seeded anywhere in the legacy
pipeline, so each fit gives a geometrically different (if topologically
similar) embedding. The actual ground-truth embeddings for Figure 3a are
ported from the legacy `reducer.pkl` files via
`scripts/port/port_umap_reducer.py` (run from that script's own pinned pixi
env, NOT this repo's `.venv`):

    cd scripts/port
    # "excl_markers" config
    pixi run python port_umap_reducer.py \\
      "/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/2-umaps/1-cafs/n_neighbors=50-min_dist=0.1-engine=umap-learn-excl_markers=dna1_dna2_fap_icsk1_icsk2_icsk3/reducer.pkl" \\
      ../../data/figures/figure3_caf_umap/excl_markers/umap_embeddings.parquet
    # "caf_markers_only" config
    pixi run python port_umap_reducer.py \\
      "/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/2-umaps/1-cafs/n_neighbors=50-min_dist=0.1-engine=umap-learn-excl_markers=beta_catenin_c_casp3_cd11b_cd20_cd3_cd31_cd4_cd44_cd45_cd66b_cd68_cd8a_dna1_dna2_e_cadherin/reducer.pkl" \\
      ../../data/figures/figure3_caf_umap/caf_markers_only/umap_embeddings.parquet

(params match this script's own N_NEIGHBORS/MIN_DIST/CONFIGS, confirmed
against `archive/scripts/02-umaps/0-umaps-cafs.py`'s `params` list -- the
legacy dir names truncate the marker list to 150 chars, hence the shorter
`caf_markers_only` name above despite excluding more markers.)
"""

from pathlib import Path

import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import NON_MARKER_CHANNELS, normalize, resolve_export_dir, resolve_output_figures_dir

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
    save_dir = resolve_output_figures_dir() / "figure3"
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
