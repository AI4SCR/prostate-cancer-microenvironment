# %%
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from ai4bmr_core.utils.logging import get_logger
from ai4bmr_core.utils.plotting import legend_from_dict
from jsonargparse import CLI

from prostate_cancer.utils import prepare_data


def main(base_dir: Path | None = None):
    logger = get_logger("phenotyping", verbose=1)

    # base_dir = base_dir or Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa')
    base_dir = base_dir or Path(f"~/data/datasets/PCa/")
    base_dir = Path(base_dir).expanduser()

    save_dir = base_dir / "02.0_clustering" / "main"
    save_dir.mkdir(parents=True, exist_ok=True)

    # FILTER
    # 01_immune
    immune = pd.read_parquet(
        base_dir / "02.0_clustering" / "immune" / f"r-{1.0}" / "data.parquet"
    )
    immune = immune[~immune.index.get_level_values("group_name").isin(["undefined"])]
    immune.index = immune.index.droplevel(
        list(set(immune.index.names) - {"object_id", "sample_name", "group_name"})
    )

    # 04_stromal
    stromal = pd.read_parquet(
        base_dir / "02.0_clustering" / "stromal" / f"r-{1.25}" / "data.parquet"
    )
    filter_ = stromal.index.get_level_values("group_name").str.startswith("stromal")
    stromal = stromal[filter_]
    stromal.index = stromal.index.droplevel(
        list(set(stromal.index.names) - {"object_id", "sample_name", "group_name"})
    )

    # 05_epithelial
    epithelial = pd.read_parquet(
        base_dir
        / "02.0_clustering"
        / "epithelial-functional"
        / f"r-{0.75}"
        / "data.parquet"
    )
    filter_ = epithelial.index.get_level_values("group_name").str.contains("epithelial")
    epithelial = epithelial[filter_]
    epithelial.index = epithelial.index.droplevel(
        list(set(epithelial.index.names) - {"object_id", "sample_name", "group_name"})
    )

    # note: remove all epithelial-basal cells, those were re-clustered
    filter_ = epithelial.index.get_level_values("group_name") == "epithelial-basal"
    epithelial = epithelial[~filter_]
    epithelial_basal = pd.read_parquet(
        base_dir / "02.0_clustering" / "epithelial-basal" / f"r-{0.5}" / "data.parquet"
    )
    epithelial_basal.index = epithelial_basal.index.droplevel(
        list(
            set(epithelial_basal.index.names)
            - {"object_id", "sample_name", "group_name"}
        )
    )

    # 06_endothelial
    endothelial = pd.read_parquet(
        base_dir / "02.0_clustering" / "endothelial" / f"r-{1.0}" / "data.parquet"
    )
    filter_ = endothelial.index.get_level_values("group_name").str.startswith(
        "endothelial"
    )
    endothelial = endothelial[filter_]
    endothelial.index = endothelial.index.droplevel(
        list(set(endothelial.index.names) - {"object_id", "sample_name", "group_name"})
    )

    # combine annotations
    cells = pd.concat([immune, stromal, endothelial, epithelial_basal, epithelial])
    # cells = pd.concat([immune, stromal, endothelial, epithelial])

    index_names = set(cells.index.names) - {"group_name"}
    annotations = (
        cells.index.to_frame().reset_index(drop=True).set_index(list(index_names))
    )

    filter_ = ~annotations.index.duplicated(keep=False)
    annotations = annotations[filter_]

    # 07_undefined
    undefined = pd.read_parquet(
        base_dir / "02.0_clustering" / "undefined" / f"r-{0.75}" / "data.parquet"
    )
    undefined.index = undefined.index.droplevel(
        list(set(undefined.index.names) - {"object_id", "sample_name", "group_name"})
    )
    undefined = (
        undefined.index.to_frame()
        .reset_index(drop=True)
        .set_index(list({"object_id", "sample_name"}))
    )

    annotations = pd.concat([annotations, undefined])
    # note: sample has been removed after initial 02.0_clustering, we need to remove it from previous annotations
    filter_ = annotations.index.get_level_values("sample_name") == "240223_012"
    annotations = annotations[~filter_]

    data = prepare_data(base_dir=base_dir)
    assert len(data) == len(annotations)
    data, annotations = data.align(annotations, join="outer", axis=0)

    parent_group_name = annotations.group_name.str.split("-").str[0]
    annotations = annotations.assign(parent_group_name=parent_group_name)
    assert annotations.isna().any().any() == False
    annotations.index = annotations.index.reorder_levels(
        ["sample_name", "object_id", "donor_block_id", "slide_code", "pat_id"]
    )

    # %%
    fname = save_dir.parent / "annotations.parquet"
    annotations.to_parquet(fname)
    logger.info(f"Saved annotations to {fname}")

    # %%
    annotations.index = annotations.index.reorder_levels(data.index.names)
    data = pd.concat([data, annotations], axis=1)
    assert data.isna().any().any() == False
    data = data.set_index(annotations.columns.tolist(), append=True)

    columns = [
        # immune
        "CD45",
        "CD20",
        "CD3",
        "CD4",
        "CD68",
        "CD8a",
        "FoxP3",
        "CD11b",
        "CD66b",
        # stromal
        "Smooth_muscle_actin",
        "Collagen1",
        "Vimentin",
        "FAP",
        "CD146",  # MCAM
        "CNN1",
        "CES1",
        "CD105",  # ENG
        "EGR1",
        "PDPN",
        # endothelial
        "CD31",
        "ERG",
        # epithelial
        "Pan-keratin",
        "E-cadherin",
        "Keratin8-18",
        "AR",
        "Synaptophysin",
        "Prostate_specific_antigen",
        "CD44",
        "Keratin5",  # basal
        "p63",  # basal
        # functional
        "Ki-67",
        "cCasp3",
        "p53",
        "YAP1",
        "beta-Catenin",
    ]

    data = data[columns]

    # %% create color_maps
    from prostate_cancer.utils import create_color_maps

    color_maps = create_color_maps(data)

    # %%
    # NOTE: we need to convert to str because of the categorical `membership` index
    row_colors = data.index.to_frame().drop(columns="object_id").astype(str)

    for label_name in row_colors.columns:
        row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

    # %%
    sort_levels = ["parent_group_name", "group_name"]
    data = data.sort_index(level=sort_levels)

    def plot_groups(data):
        cg = sns.clustermap(
            data,
            row_colors=row_colors,
            row_cluster=False,
            col_cluster=False,
            figsize=(15, 10),
        )

        cg.ax_heatmap.set_yticklabels([])
        cg.ax_heatmap.set_ylabel("objects")

        legend_elements = legend_from_dict(color_maps["group_name"])
        cg.ax_heatmap.legend(
            handles=legend_elements,
            loc="upper left",
            bbox_to_anchor=(1.05, 1),
            fontsize=6,
        )

        legend_elements = legend_from_dict(color_maps["parent_group_name"])
        cg.ax_row_colors.legend(
            handles=legend_elements,
            loc="lower right",
            bbox_to_anchor=(0, 0),
            fontsize=6,
        )

        cg.figure.tight_layout()
        # cg.figure.show()
        return cg

    cg = plot_groups(data)
    cg.figure.savefig(save_dir / "groups.png", dpi=300)
    plt.close(cg.figure)

    # %% group medians
    def plot_medians(data, level, row_cluster: bool = False):
        pdat = data.groupby(level, observed=True).median()
        pdat = pdat.sort_values(level)

        row_colors = pdat.index.to_frame().astype(str)

        for label_name in row_colors.columns:
            row_colors[label_name] = row_colors[label_name].map(color_maps[label_name])

        cg = sns.clustermap(
            pdat,
            row_colors=row_colors,
            row_cluster=row_cluster,
            col_cluster=False,
            figsize=(15, 10),
        )

        cg.ax_heatmap.set_yticklabels([])
        cg.ax_heatmap.set_ylabel("objects")

        legend_elements = legend_from_dict(color_maps["group_name"])
        cg.ax_heatmap.legend(
            handles=legend_elements,
            loc="upper left",
            bbox_to_anchor=(1.05, 1),
            fontsize=6,
        )

        legend_elements = legend_from_dict(color_maps["parent_group_name"])
        cg.ax_row_colors.legend(
            handles=legend_elements,
            loc="lower right",
            bbox_to_anchor=(0, 0),
            fontsize=6,
        )

        cg.figure.tight_layout()
        return cg

    cg = plot_medians(data, ["parent_group_name", "group_name"], row_cluster=True)
    cg.figure.savefig(save_dir / "groups-medians.png", dpi=300)
    plt.close(cg.figure)

    # %% group medians parent
    cg = plot_medians(data, "parent_group_name", row_cluster=True)
    cg.figure.savefig(save_dir / "parent-groups-medians.png", dpi=300)
    plt.close(cg.figure)

    # %%
    for grp_name, grp_data in data.groupby("parent_group_name", observed=True):
        cg = plot_groups(grp_data)
        cg.figure.savefig(save_dir / f"{grp_name}.png", dpi=300)
        plt.close(cg.figure)

        cg = plot_medians(
            grp_data, ["parent_group_name", "group_name"], row_cluster=True
        )
        cg.figure.savefig(save_dir / f"{grp_name}-medians.png", dpi=300)
        plt.close(cg.figure)


if __name__ == "__main__":
    CLI(main)
