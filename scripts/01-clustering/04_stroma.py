# %%
import pickle
from pathlib import Path

import pandas as pd
import seaborn as sns
from loguru import logger as base_logger
from jsonargparse import CLI
from matplotlib import pyplot as plt
from prostate_cancer.plotting import legend_from_dict

from prostate_cancer.cluster import cluster
from prostate_cancer.utils import prepare_data


def main(base_dir: Path | None = None, resolution: float = 1.25):
    logger = base_logger.bind(task="phenotyping")

    markers = [
        # epithelial
        "Pan-keratin",
        "E-cadherin",
        "Keratin8-18",
        # endothelial
        "CD31",
        # immune
        "CD45",
        # CAFs
        "Smooth_muscle_actin",
        "Collagen1",
        "Vimentin",
        "CD146",  # MCAM
        "CNN1",
        "CES1",
        "CD105",  # ENG
        "EGR1",
        "PDPN",
        "AR",
    ]
    marker_groups = {
        "epithelial": ["Pan-keratin", "E-cadherin", "Keratin8-18"],
        "endothelial": ["CD31", "CD146", "CD105"],
        "stromal": ["Vimentin", "Smooth_muscle_actin", "Collagen1"],
    }

    base_dir = base_dir or Path("/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa")
    # base_dir = base_dir or Path('~/data/datasets/PCa')
    # base_dir = base_dir or Path('~/data/pca-v3')
    base_dir = Path(base_dir).expanduser()

    save_dir = base_dir / "02.0_clustering" / "stromal" / f"r-{resolution}"
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
        cells = pd.read_parquet(
            base_dir
            / "02.0_clustering"
            / "epithelial-non-epithelial"
            / "r-1.0"
            / "data.parquet"
        )
        cells = cells[
            cells.index.get_level_values("group_name").str.contains("stromal")
        ]

        endothelial = pd.read_parquet(
            base_dir
            / "02.0_clustering"
            / "immune-non-immune"
            / "r-1.5"
            / "data.parquet"
        )
        endothelial = endothelial[
            endothelial.index.get_level_values("group_name") == "endothelial"
        ]

        # NOTE: this will introduce duplicated cells as we have the undefined group included as well when 02.0_clustering
        #  the epithelial-non-epithelial data
        undefined = pd.read_parquet(
            base_dir / "02.0_clustering" / "immune" / "r-1.0" / "data.parquet"
        )
        undefined = undefined[
            undefined.index.get_level_values("group_name").isin(["undefined"])
        ]

        cells = pd.concat((cells, undefined, endothelial))

        redundant_index_levels = set(cells.index.names) - set(data.index.names)
        include_index = cells.index.droplevel(list(redundant_index_levels))
        assert include_index.duplicated().any() == True
        include_index = include_index.drop_duplicates()
        data = data.loc[include_index]

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

    # from ai4bmr_core.utils.sampling import sample_min_per_group_then_uniform
    # grouped = data.groupby('membership', observed=True)
    # pdat = sample_min_per_group_then_uniform(grouped, n=int(1e5))
    # pdat = pdat.sample(frac=1)  # shuffle for plotting
    # assert pdat.index.duplicated().any() == False
    #
    # embedding = pd.DataFrame(embedding, index=data.index, columns=['umap1', 'umap2'])
    # embedding = embedding.loc[pdat.index].values

    # %%
    # from 02.0_clustering.utils import plot_umap_index
    # plot_umap_index(pdat, embedding=embedding, color_maps=color_maps, save_dir=save_dir)

    # %%
    # from 02.0_clustering.utils import plot_umap_columns
    # plot_umap_columns(pdat, embedding=embedding, save_dir=save_dir)

    # %% PRE-FILTER FOR VIMENTIN
    # save_dir = base_dir / '02.0_clustering' / 'caf-non-caf-subset' / f'r-{resolution}'
    # save_dir.mkdir(exist_ok=True, parents=True)
    # data = data[data.index.get_level_values('membership').isin(['3', '4', '5', '6', '8'])]

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

    for col in set(pdat) - {"proportion"}:
        ax = pdat[[col, "proportion"]].sort_values(col).plot(kind="bar")
        ax.figure.savefig(save_dir / f"median_{col}.png", dpi=300)
        plt.close(ax.figure)

    # %% MEDIAN EXPRESSION OF SETS
    for name, col in marker_groups.items():
        ax = pdat[col + ["proportion"]].sort_values(col).plot(kind="bar")
        ax.figure.savefig(save_dir / f"median_{name}.png", dpi=300)
        plt.close(ax.figure)

    # %%
    ax = pdat.plot(kind="bar")
    # ax.figure.show()
    ax.figure.savefig(save_dir / "median_markers.png", dpi=300)
    plt.close(ax.figure)

    # %% create color_maps
    from prostate_cancer.utils import create_color_maps

    color_maps = create_color_maps(data)

    # %% CLUSTERMAP MEDIAN's
    row_colors = median_.index.to_frame()

    for label_name in row_colors.columns:
        row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

    median_ = median_.sort_index(level="membership")
    cg = sns.clustermap(
        median_, row_colors=row_colors, row_cluster=True, col_cluster=False
    )

    cg.ax_heatmap.set_yticklabels([])
    cg.ax_heatmap.set_ylabel("objects")

    cluster_ids = (
        data.index.get_level_values("membership")
        .unique()
        .astype(int)
        .sort_values()
        .astype(str)
    )
    legend_dict = {i: color_maps["membership"][i] for i in cluster_ids}
    legend_elements = legend_from_dict(legend_dict)
    cg.ax_row_colors.legend(
        handles=legend_elements, loc="lower right", bbox_to_anchor=(0, 0), fontsize=6
    )

    # cg.figure.show()
    cg.figure.savefig(save_dir / "clustermap-median.png", dpi=300)
    plt.close(cg.figure)

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

    cluster_ids = (
        data.index.get_level_values("membership")
        .unique()
        .astype(int)
        .sort_values()
        .astype(str)
    )
    legend_dict = {i: color_maps["membership"][i] for i in cluster_ids}
    legend_elements = legend_from_dict(legend_dict)
    cg.ax_row_colors.legend(
        handles=legend_elements, loc="lower right", bbox_to_anchor=(0, 0), fontsize=6
    )

    # cg.figure.show()
    cg.figure.tight_layout()
    cg.figure.savefig(save_dir / "clustermap.png", dpi=300)
    plt.close(cg.figure)

    # %%
    logger.info("completed")


if __name__ == "__main__":
    CLI(main)
