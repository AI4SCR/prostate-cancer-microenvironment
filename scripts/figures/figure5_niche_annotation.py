# %%
"""Reproduce Figure 5's niche naming/annotation step.

1:1 port of the old repo's `000_paper/11_niches/110_analysis/01_annotation_v2.py`
(paths only changed; see `figure_script_mapping.md`) -- maps each of the 24
raw k-means cluster IDs (`figure5_niche_clustering.py`'s output) to a
human-readable niche name and a meta-niche grouping, via a manually-curated
lookup table.

BLOCKED: that lookup table, `niche_annotations_revised.xlsx`, is genuinely
missing -- not just permission-denied like the rest of this figure's
dependencies. Confirmed absent from every location checked (the live but
inaccessible `/users/mensmeng/...` path, the consolidated `LEGACY_DATA_DIR`
copy, the pre-migration repo, and the shared `/work/.../prometex/data/PCa/`
tree). See `missing_files.md`. This script will raise `FileNotFoundError`
at the `pd.read_excel()` call below until that file is supplied.

Its OWN output, however, was already computed once (presumably by whoever
had access to the xlsx) and is staged, precomputed, at
`$LEGACY_DATA_DIR/5-niches/annotation/{niche_annotations_v2.csv,
clusters_annotated_v2.parquet}` -- downstream figure scripts (heatmap,
correlation, Figure 6/7) read that precomputed file directly rather than
depending on this script completing.
"""
from pathlib import Path

import pandas as pd
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import resolve_export_dir, resolve_legacy_dir


def main(export_dir: Path | None = None, legacy_dir: Path | None = None):
    export_dir = export_dir or resolve_export_dir()
    legacy_dir = legacy_dir or resolve_legacy_dir()
    save_dir = export_dir / "figures" / "figure5"
    save_dir.mkdir(parents=True, exist_ok=True)

    # %% figure5_niche_clustering.py's output
    df_clusters = pd.read_parquet(save_dir / "clusters.parquet", engine="fastparquet")
    logger.info(f"clusters shape: {df_clusters.shape}")

    # %% GENUINELY MISSING -- see module docstring and missing_files.md
    annot_excel = legacy_dir / "PCA_NHOODs_clean" / "niche_annotations_revised.xlsx"
    df_anno = pd.read_excel(annot_excel)

    annotation_dict_niche = dict(zip(df_anno["cluster"].astype(str), df_anno["niche"]))
    annotation_dict_meta_niche = dict(zip(df_anno["niche"], df_anno["meta_niche"]))
    color_dict_niche = dict(zip(df_anno["niche"], df_anno["niche_color"]))
    color_dict_meta_niche = dict(zip(df_anno["meta_niche"], df_anno["meta_niche_color"]))

    cluster_name = df_clusters.columns[0]
    df_clusters[cluster_name] = df_clusters[cluster_name].astype(str)
    df_clusters["niche"] = df_clusters[cluster_name].map(annotation_dict_niche).fillna("unassigned")

    # %%
    df_annotation_niche = pd.DataFrame.from_dict(annotation_dict_niche, orient="index", columns=["niche"])
    df_annotation_niche.index.name = "cluster"
    df_annotation_niche.reset_index(inplace=True)
    df_annotation_niche["meta_niche"] = df_annotation_niche["niche"].map(annotation_dict_meta_niche).fillna("unassigned")
    df_annotation_niche["niche_color"] = df_annotation_niche["niche"].map(color_dict_niche).fillna("#7f7f7f")
    df_annotation_niche["meta_niche_color"] = df_annotation_niche["meta_niche"].map(color_dict_meta_niche).fillna("#7f7f7f")
    df_annotation_niche.to_csv(save_dir / "niche_annotations_v2.csv", index=False)

    # %% legacy repeats this exact assignment a second time (its lines 134-136,
    # byte-identical to the one above) before the final save -- reproduced
    # here even though idempotent, for strict fidelity to the original script.
    cluster_name = df_clusters.columns[0]
    df_clusters[cluster_name] = df_clusters[cluster_name].astype(str)
    df_clusters["niche"] = df_clusters[cluster_name].map(annotation_dict_niche).fillna("unassigned")

    df_clusters["meta_niche"] = df_clusters["niche"].map(annotation_dict_meta_niche).fillna("unassigned")
    df_clusters["niche"] = df_clusters["niche"].astype("category")
    df_clusters["meta_niche"] = df_clusters["meta_niche"].astype("category")
    df_clusters.to_parquet(save_dir / "clusters_annotated_v2.parquet", engine="fastparquet", index=True)
    logger.info(f"saved annotated clusters to {save_dir}")


if __name__ == "__main__":
    CLI(main)
