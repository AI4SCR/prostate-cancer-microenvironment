# %%
import pickle
from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns
from ai4bmr_core.utils.logging import get_logger
from jsonargparse import CLI

from prostate_cancer.cluster import cluster
from prostate_cancer.utils import prepare_data


def main(base_dir: Path | None = None, resolution: float = 1.5):
    logger = get_logger("phenotyping", verbose=1)

    markers = [
        "CD20",
        "CD3",
        "CD4",
        "CD68",
        "CD8a",
        "FoxP3",
        "CD45",
        "CD11b",
        "CD66b",  # immune
        "CD31",  # endothelial
        "E-cadherin",
        "Pan-keratin",  # epithelial
        # 'Smooth_muscle_actin', 'Vimentin', 'Collagen1'  # fibroblast
        # 'Collagen1', 'Smooth_muscle_actin' # fibroblast
    ]
    marker_groups = (
        ["Pan-keratin", "E-cadherin"],
        ["CD45", "CD31", "Pan-keratin", "E-cadherin"],
        ["CD45", "FoxP3"],
    )

    base_dir = base_dir or Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa")
    base_dir = Path(base_dir).expanduser()

    save_dir = base_dir / "02.0_clustering" / "immune-non-immune" / f"r-{resolution}"
    save_dir.mkdir(exist_ok=True, parents=True)

    # %%
    data_path = save_dir / "data.pkl"

    if data_path.exists():
        logger.info(f"loading data from {data_path}")
        with open(data_path, "rb") as f:
            data = pickle.load(f)
            data, embedding, result = data["data"], data["embedding"], data["result"]
    else:
        logger.info(f"02.0_clustering")
        data = prepare_data(base_dir=base_dir)

        # FILTER
        data = data[markers]

        # CLUSTER
        data, embedding, result = cluster(
            data,
            resolution=resolution,
            graph_method="scanpy",
            cluster_method="igraph",
            n_neighbors=15,
        )
        if data.index.get_level_values("membership").nunique() > 50:
            raise ValueError("too many clusters")

        with open(data_path, "wb") as f:
            d = dict(data=data, embedding=embedding, result=result)
            pickle.dump(d, f)

    # %%
    logger.info("creating figures")

    # %% MEDIAN EXPRESSION
    median_ = data.groupby("membership", observed=True).median()
    median_.index = median_.index.astype(str)

    pdat = pd.concat(
        (
            median_,
            data.index.get_level_values("membership").value_counts(normalize=True),
        ),
        axis=1,
    )

    median_cd45 = median_["CD45"]
    cd45_threshold = 0.25
    for col in set(pdat) - {"proportion"}:
        ax = pdat[[col, "proportion"]].sort_values(col).plot(kind="bar")

        if col == "CD45":
            ax.axhline(cd45_threshold, color="red", linestyle="--")

        ax.figure.savefig(save_dir / f"median_{col}.png", dpi=300)

    # %% MEDIAN EXPRESSION OF SETS
    for col in marker_groups:
        ax = pdat[col + ["proportion"]].sort_values(col).plot(kind="bar")
        name = "_".join(col)
        ax.figure.savefig(save_dir / f"median_{name}.png", dpi=300)

    # %%
    ax = pdat.plot(kind="bar")
    # ax.figure.show()
    ax.figure.savefig(save_dir / "median_markers.png", dpi=300)

    # %% create color_maps
    from prostate_cancer.utils import create_color_maps

    color_maps = create_color_maps(data)

    # %% CLUSTERMAP MEDIAN's
    row_colors = median_.index.to_frame()

    for label_name in row_colors.columns:
        row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

    median_ = median_.sort_index(level="membership")
    median_ = median_.sort_values(["CD45", "FoxP3"])
    cg = sns.clustermap(median_, row_colors=row_colors, row_cluster=False)

    cg.ax_heatmap.set_yticklabels([])
    cg.ax_heatmap.set_ylabel("objects")
    # cg.figure.show()
    cg.figure.savefig(save_dir / "clustermap-median.png", dpi=300)

    # %%
    # NOTE: we need to convert to str because of the categorical `membership` index
    row_colors = data.index.to_frame().drop(columns="object_id").astype(str)

    # row_colors.replace(color_maps) # does not work
    for label_name in row_colors.columns:
        row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

    # %%
    cg = sns.clustermap(
        data.sort_index(level="membership"),
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=False,
    )

    cg.ax_heatmap.set_yticklabels([])
    cg.ax_heatmap.set_ylabel("objects")
    # cg.figure.show()
    cg.figure.tight_layout()
    cg.figure.savefig(save_dir / "clustermap.png", dpi=300)

    # %%
    cmap = matplotlib.colormaps.get_cmap("plasma")

    row_colors["median_cd45"] = median_cd45.loc[
        row_colors.index.get_level_values("membership")
    ].values
    row_colors["is_immune"] = (row_colors["median_cd45"] >= cd45_threshold).map(
        {True: "green", False: "red"}
    )
    row_colors["median_cd45"] = row_colors["median_cd45"].map(cmap).map(list)

    median_pan_keratin = median_["Pan-keratin"]
    row_colors["median_pan_keratin"] = median_pan_keratin.loc[
        row_colors.index.get_level_values("membership")
    ].values
    row_colors["median_pan_keratin"] = (
        row_colors["median_pan_keratin"].map(cmap).map(list)
    )

    median_cd31 = median_["CD31"]
    row_colors["median_cd31"] = median_cd31.loc[
        row_colors.index.get_level_values("membership")
    ].values
    row_colors["median_cd31"] = row_colors["median_cd31"].map(cmap).map(list)

    # order memberships by median CD45
    membership = pd.Categorical(
        data.index.get_level_values("membership"),
        categories=median_cd45.sort_values().index,
        ordered=True,
    )
    pdat = (
        data.droplevel("membership")
        .assign(membership=membership)
        .set_index("membership", append=True)
        .sort_index(level="membership")
    )

    cg = sns.clustermap(
        pdat.sort_index(level="membership"),
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=False,
    )

    cg.ax_heatmap.set_yticklabels([])
    cg.ax_heatmap.set_ylabel("objects")
    # cg.figure.show()
    cg.figure.tight_layout()
    cg.figure.savefig(save_dir / "clustermap-medians.png", dpi=300)

    # %%
    logger.info("completed")


if __name__ == "__main__":
    CLI(main)
