# %%
import pickle
from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns
from ai4bmr_core.utils.logging import get_logger
from ai4bmr_core.utils.plotting import legend_from_dict
from jsonargparse import CLI


def main(base_dir: Path | None = None, resolution: float = 1.5):
    logger = get_logger("phenotyping", verbose=1)

    # base_dir = base_dir or Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa')
    base_dir = base_dir or Path("~/data/datasets/PCa")
    # base_dir = base_dir or Path('~/data/pca-v3')
    base_dir = Path(base_dir).expanduser()

    save_dir = base_dir / "02.0_clustering" / "immune-non-immune" / f"r-{resolution}"
    assert save_dir.exists()

    # %%
    data_path = save_dir / "data.pkl"
    logger.info(f"loading data from {data_path}")
    with open(data_path, "rb") as f:
        data = pickle.load(f)
        data, embedding, result = data["data"], data["embedding"], data["result"]

    # %%
    cluster_id_2_name = {
        "18": "immune",
        "19": "immune",
        "20": "immune",
        "10": "immune",
        "2": "immune",
        "8": "immune",
        "1": "immune",
        "9": "immune-epithelial",
        "0": "immune",
        "12": "immune",  # low cd45 signal
        "25": "endothelial",
        "5": "endothelial",
        "39": "low-signal",
        "32": "low-signal",
        "40": "low-signal",
        "21": "low-signal",
        "22": "low-signal",
        "23": "low-signal",
        "34": "low-signal",
        "36": "low-signal",
        "37": "low-signal",
        "33": "low-signal",
        "38": "low-signal",
        "41": "low-signal",
        "31": "low-signal",
        "16": "low-signal",
        "35": "low-signal",
        "26": "low-signal",
        "13": "low-signal",
        "11": "low-signal",
        "4": "low-signal",
        "29": "epithelial",
        "27": "epithelial",
        "28": "epithelial",
        "30": "epithelial",
        "17": "epithelial",
        "14": "epithelial",
        "6": "epithelial",
        "3": "epithelial",
        "15": "epithelial",
        "7": "epithelial",
        "24": "epithelial",
    }
    assert set(cluster_id_2_name) == set(
        data.index.get_level_values("membership").unique()
    )

    cmap = {
        "low-signal": matplotlib.colors.to_rgb("black"),
        "immune": matplotlib.colors.to_rgb("darkgreen"),
        "immune-epithelial": matplotlib.colors.to_rgb("green"),
        "endothelial": matplotlib.colors.to_rgb("brown"),
        "epithelial": matplotlib.colors.to_rgb("blue"),
    }

    # %%
    group_name = pd.Categorical(
        data.index.get_level_values("membership").map(lambda x: cluster_id_2_name[x]),
    )
    assert group_name.isna().sum() == 0

    data = (
        data.assign(group_name=group_name)
        .set_index("group_name", append=True)
        .sort_index(level="group_name")
    )

    # %% create color_maps
    from prostate_cancer.utils import create_color_maps

    color_maps = create_color_maps(data)
    color_maps["group_name"] = cmap

    # %%
    # NOTE: we need to convert to str because of the categorical `membership` index
    row_colors = data.index.to_frame().drop(columns="object_id").astype(str)

    row_colors["is_immune"] = row_colors["group_name"].str.contains("immune")
    color_maps["is_immune"] = {True: "green", False: "red"}

    for label_name in row_colors.columns:
        row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

    # %%
    plasma = matplotlib.colormaps.get_cmap("plasma")

    median_ = data.groupby("membership", observed=True).median()
    median_.index = median_.index.astype(str)

    median_cd45 = median_["CD45"]
    row_colors["median_cd45"] = median_cd45.loc[
        row_colors.index.get_level_values("membership")
    ].values
    row_colors["median_cd45"] = row_colors["median_cd45"].map(plasma).map(list)

    median_pan_keratin = median_["Pan-keratin"]
    row_colors["median_pan_keratin"] = median_pan_keratin.loc[
        row_colors.index.get_level_values("membership")
    ].values
    row_colors["median_pan_keratin"] = (
        row_colors["median_pan_keratin"].map(plasma).map(list)
    )

    median_cd31 = median_["CD31"]
    row_colors["median_cd31"] = median_cd31.loc[
        row_colors.index.get_level_values("membership")
    ].values
    row_colors["median_cd31"] = row_colors["median_cd31"].map(plasma).map(list)

    # %%
    # NOTE: create categories ordered by the median CD45 expression
    membership = pd.Categorical(
        data.index.get_level_values("membership"),
        categories=median_cd45.sort_values().index,
        ordered=True,
    )

    data = data.assign(m=membership).sort_values("m").drop(columns="m")

    # %%
    cg = sns.clustermap(
        data, row_colors=row_colors, row_cluster=False, col_cluster=False
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

    # %%
    data.to_parquet(save_dir / "data.parquet")
    logger.info("completed")


if __name__ == "__main__":
    CLI(main)
