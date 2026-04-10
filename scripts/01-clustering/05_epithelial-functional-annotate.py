# %%
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from ai4bmr_core.utils.logging import get_logger
from ai4bmr_core.utils.plotting import legend_from_dict
from jsonargparse import CLI


def main(base_dir: Path | None = None):
    logger = get_logger("phenotyping", verbose=1)

    # base_dir = base_dir or Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa')
    base_dir = base_dir or Path(f"~/data/datasets/PCa/")
    # base_dir = base_dir or Path('~/data/pca-v3')
    base_dir = Path(base_dir).expanduser()

    save_dir = base_dir / "02.0_clustering" / "epithelial-functional" / f"r-{0.75}"
    assert save_dir.exists()

    # %%
    data_path = save_dir / "data.pkl"

    logger.info(f"loading data from {data_path}")
    with open(data_path, "rb") as f:
        data = pickle.load(f)
        data, embedding, result = data["data"], data["embedding"], data["result"]

    # %%
    index = data.index.to_frame().reset_index(drop=True)
    index[index.sample_name == "240223_012"].membership.value_counts()
    index[index.membership == "17"].sample_name.value_counts()
    index.groupby("membership", observed=True).size().loc["17"]

    # %%
    cluster_id_2_name = {
        "0": "epithelial-basal",  # basal?
        "1": "undefined",  # expresses PSA which can be luminal but also p63 which is basal, everything else is low
        "2": "epithelial-Synaptophysin_pos",
        "3": "epithelial-ERG_pos",
        "4": "endothelial-Vimentin_neg",
        "5": "epithelial",
        "6": "epithelial-PSA_pos",
        "7": "epithelial-PSA_pos",
        "8": "epithelial",
        "9": "epithelial-ERG_pos",
        "10": "epithelial",
        "11": "epithelial",
        "12": "epithelial-ERG_pos-CD44_pos",
        "13": "epithelial",
        "14": "epithelial-ERG_pos-p53_pos",
        "15": "epithelial-ERG_pos",
        "16": "undefined",
        "17": "epithelial-Ki67_pos",
        "18": "endothelial-ERG_pos",
        "19": "endothelial",
        "20": "epithelial-ERG_pos-PSA_pos",
        "21": "epithelial-p53_pos",
        "22": "epithelial-PSA_pos-ERG_pos",
        "23": "epithelial-Synaptophysin_pos-p53_pos",  # collagen1_high but probably contaminated
        "24": "epithelial-ERG_pos-p53_pos",
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
