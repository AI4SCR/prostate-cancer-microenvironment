# %%
import pickle
from pathlib import Path

import pandas as pd
import seaborn as sns
from ai4bmr_core.utils.logging import get_logger
from ai4bmr_core.utils.plotting import legend_from_dict
from jsonargparse import CLI

from prostate_cancer.utils import prepare_data


def main(base_dir: Path | None = None):
    logger = get_logger("phenotyping", verbose=1)

    # base_dir = base_dir or Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa')
    base_dir = base_dir or Path("~/data/datasets/PCa")
    # base_dir = base_dir or Path('~/data/pca-v3')
    base_dir = Path(base_dir).expanduser()

    save_dir = base_dir / "02.0_clustering" / "immune" / f"r-{1.0}"
    assert save_dir.exists()

    # %%
    data_path = save_dir / "data.pkl"

    logger.info(f"loading data from {data_path}")
    with open(data_path, "rb") as f:
        data = pickle.load(f)
        data, embedding, result = data["data"], data["embedding"], data["result"]
        markers = data.columns.tolist()

    membership = data.reset_index("membership")["membership"]

    # %%
    base_dir = base_dir or Path("~/data/datasets/PCa")
    base_dir = Path(base_dir).expanduser()
    data = prepare_data(base_dir)

    data = (
        pd.concat((data, membership), axis=1)
        .dropna()
        .set_index("membership", append=True)
    )
    # data = data.loc[:, ~data.columns.isin(['DNA1', 'DNA2', 'ICSK1', 'ICSK2', 'ICSK3'])]
    data = data[markers + ["CD45"]]

    # %%
    # NOTE: should we label the t_cells_helper_missing_cd3? Would be cluster 15, 0, 1
    cluster_id_2_name = {
        "23": "immune-t_cells_helper",
        "17": "immune-t_cells_helper_cytotoxic_mix",
        "5": "immune-t_cells_cytotoxic",
        "18": "immune-t_cells",
        "12": "immune-b_cells_t_cells_helper_mix",
        "22": "immune-t_cells_regulatory",
        "13": "immune-t_cells_helper",
        "4": "immune-t_cells_helper",
        "11": "immune-t_cells_helper_macrophages_mix",
        "3": "immune-t_cells_cytotoxic",
        "9": "immune-t_cells_helper",
        "19": "immune-b_cells",
        "2": "immune-macrophages",
        "0": "immune-t_cells_cytotoxic",
        "7": "immune-t_cells_helper",
        "20": "undefined",  # b_cells -> undefined
        "10": "immune-macrophages",
        "6": "undefined",  # NOTE: t_cells_helper ? large fraction, relative expression indicates t_helper but no CD3 signal
        "14": "immune-macrophages",
        "21": "immune-neutrophils",
        "15": "immune-macrophages",
        "1": "undefined",  # t_cells_cytotoxic ?
        "16": "immune-neutrophils",
        "8": "immune-macrophages",  # NOTE: low signal but 4th for C11b and 8th for CD68
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
        data,
        figsize=(20, 10),
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

    # %%
    data.to_parquet(save_dir / "data.parquet")
    logger.info("completed")


if __name__ == "__main__":
    CLI(main)
