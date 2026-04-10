# %%
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from loguru import logger as base_logger
from jsonargparse import CLI
from prostate_cancer.plotting import legend_from_dict


def main(base_dir: Path | None = None):
    logger = base_logger.bind(task="phenotyping")

    # base_dir = base_dir or Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa')
    base_dir = base_dir or Path(f"~/data/datasets/PCa/")
    # base_dir = base_dir or Path('~/data/pca-v3')
    base_dir = Path(base_dir).expanduser()

    save_dir = base_dir / "02.0_clustering" / "epithelial-basal" / f"r-{0.5}"
    assert save_dir.exists()

    # %%
    data_path = save_dir / "data.pkl"

    logger.info(f"loading data from {data_path}")
    with open(data_path, "rb") as f:
        data = pickle.load(f)
        data, embedding, result = data["data"], data["embedding"], data["result"]

    # %%
    cluster_id_2_name = {
        "0": "epithelial-luminal",
        "1": "epithelial-luminal",
        "2": "epithelial-luminal",
        "3": "epithelial-luminal",
        "4": "epithelial-luminal",
        "5": "epithelial-basal",
        "6": "epithelial-transient",
        "7": "epithelial-basal",  # AR-, KRT8-18- but borderline p63 and low KRT5, in general low but mixed signal
        "8": "epithelial-basal",
        "9": "epithelial-basal",
        "10": "epithelial-basal",
        "11": "epithelial-basal",
        "12": "epithelial-basal",
        "13": "epithelial-basal",
        "14": "epithelial-basal",
        "15": "epithelial-basal",
        "16": "epithelial-basal",
    }

    assert set(cluster_id_2_name) == set(
        data.index.get_level_values("membership").unique()
    )

    import colorcet

    cmap = {
        k: colorcet.glasbey_light[i]
        for i, k in enumerate(sorted(set(cluster_id_2_name.values())))
    }
    cmap["undefined"] = "black"

    # %%
    group_name = pd.Categorical(
        data.index.get_level_values("membership").map(lambda x: cluster_id_2_name[x]),
    )
    assert group_name.isna().sum() == 0

    data = (
        data.assign(group_name=group_name)
        .set_index("group_name", append=True)
        .sort_index(level=["group_name", "membership"])
    )

    # %% create color_maps
    from prostate_cancer.utils import create_color_maps

    color_maps = create_color_maps(data)
    color_maps["group_name"] = cmap

    # %%
    # NOTE: we need to convert to str because of the categorical `membership` index
    row_colors = data.index.to_frame().drop(columns="object_id").astype(str)
    color_maps["group_name"] = cmap

    for label_name in row_colors.columns:
        row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

    # %%
    cg = sns.clustermap(
        data.sort_values("group_name"),
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=False,
    )

    cg.ax_heatmap.set_yticklabels([])
    cg.ax_heatmap.set_ylabel("objects")

    legend_elements = legend_from_dict(cmap)
    cg.ax_heatmap.legend(
        handles=legend_elements, loc="upper left", bbox_to_anchor=(1.05, 1), fontsize=6
    )

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

    cg.figure.tight_layout()
    # cg.figure.show()
    cg.figure.savefig(save_dir / "clustermap-groups.png", dpi=300)
    plt.close(cg.figure)

    # %% group medians
    pdat = data.groupby("group_name", observed=True).median()
    row_colors = pdat.index.to_frame().astype(str)
    color_maps["group_name"] = cmap

    for label_name in row_colors.columns:
        row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

    cg = sns.clustermap(
        pdat.sort_values("group_name"),
        row_colors=row_colors,
        row_cluster=True,
        col_cluster=True,
    )

    cg.ax_heatmap.set_yticklabels([])
    cg.ax_heatmap.set_ylabel("objects")

    legend_elements = legend_from_dict(cmap)
    cg.ax_heatmap.legend(
        handles=legend_elements, loc="upper left", bbox_to_anchor=(1.05, 1), fontsize=6
    )

    cg.figure.tight_layout()
    # cg.figure.show()
    cg.figure.savefig(save_dir / "clustermap-groups-medians.png", dpi=300)
    plt.close(cg.figure)

    # %%
    data.to_parquet(save_dir / "data.parquet")
    logger.info("completed")


if __name__ == "__main__":
    CLI(main)
