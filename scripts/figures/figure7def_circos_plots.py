from pathlib import Path
from typing import Dict

import matplotlib
import numpy as np
import pandas as pd
import yaml
from jsonargparse import CLI
from loguru import logger
from matplotlib import cm, colors
from pycirclize import Circos

from prostate_cancer.utils import resolve_data_dir, resolve_output_figures_dir

matplotlib.use("Agg")

NICHES = {
    "d": "luminal_infiltrated",  # niche 2
    "e": "tumor_CAF1_lymphocytes",  # niche 8
    "f": "luminal_CAF1(CD105High)",  # niche 9
}


def check_triangular_zero_xor(df: pd.DataFrame):
    """Checks if the strictly upper OR strictly lower triangle of a square
    DataFrame is all zeros, but not both."""
    if df.shape[0] != df.shape[1]:
        return False
    data = df.values
    is_upper_zero = np.all(np.triu(data, k=1) == 0)
    is_lower_zero = np.all(np.tril(data, k=-1) == 0)
    return is_upper_zero ^ is_lower_zero


def get_color_map(df: pd.DataFrame) -> Dict[str, str]:
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "
    if df.index.names != ["cell_type_1", "cell_type_2"]:
        df.index.set_names(["cell_type_1", "cell_type_2"], inplace=True)
    df_new = df.reset_index()
    labels = list(set(df_new["cell_type_1"].unique()).union(set(df_new["cell_type_2"].unique())))
    cmap = cm.get_cmap("tab10")
    color_indices = {label: i for i, label in enumerate(labels)}
    color_map = {ct: cmap(color_indices[ct]) for ct in labels}
    return color_map


def data_for_circos_plot(df: pd.DataFrame, color: str):
    assert color in ["aggregated_interaction", "above_median_fraction"]
    assert "above_median_fraction" and "aggregated_interaction" in df.columns
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    if df.index.names != ["cell_type_1", "cell_type_2"]:
        df.index.set_names(["cell_type_1", "cell_type_2"], inplace=True)

    # WIDTH
    width = [col for col in df.columns.tolist() if col != color][0]
    width_df = df[width]
    width_df = width_df.reset_index()
    width_df.columns = ["from", "to", "Value"]
    width_df["Value"].dropna()
    width_dict = {(row["from"], row["to"]): row["Value"] for _, row in width_df.iterrows()}

    # COLOR
    color_df = df[color]
    color_df = color_df.reset_index()
    color_df.columns = ["from", "to", "Value"]
    color_df["Value"].dropna()
    color_dict = {(row["from"], row["to"]): row["Value"] for _, row in color_df.iterrows()}

    # SECTORS
    sectors_df = df[width].copy()
    sectors_df = sectors_df.reset_index()
    sectors_df.columns = ["cell_type_1", "cell_type_2", "width"]
    sectors_df = sectors_df.pivot(index="cell_type_2", columns="cell_type_1", values="width")
    sectors_df = sectors_df.fillna(0)
    assert check_triangular_zero_xor(sectors_df), "Either the upper or lower triangle values of sectors_df must be 0, but not both."

    return {"color_dict": color_dict, "width_dict": width_dict, "sectors_df": sectors_df}


def get_circos(df: pd.DataFrame, color: str = "aggregated_interaction", color_map: Dict[str, str] | None = None):
    assert color in ["aggregated_interaction", "above_median_fraction"]
    assert "above_median_fraction" and "aggregated_interaction" in df.columns
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    dicts = data_for_circos_plot(df=df, color=color)
    color_dict = dicts["color_dict"]
    width_dict = dicts["width_dict"]
    sectors_df = dicts["sectors_df"]

    cell_color_map = get_color_map(df) if color_map is None else color_map
    color_var = color

    def link_handler(from_label, to_label):
        val = color_dict.get((from_label, to_label)) or color_dict.get((to_label, from_label), 0)
        lw = width_dict.get((from_label, to_label)) or width_dict.get((to_label, from_label), 1)

        filt = val if color_var == "aggregated_interaction" else lw

        if filt <= 0.001:
            color = "none"
            lw = 0
        else:
            cmap = cm.get_cmap("Reds")
            val_min = min(color_dict.values())
            val_max = max(color_dict.values())
            norm = colors.Normalize(vmin=val_min, vmax=val_max)
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            color = sm.to_rgba(val)

        return dict(ec="none", lw=lw, fc=color, alpha=0.7)

    circos = Circos.chord_diagram(
        sectors_df,
        space=2,
        cmap=cell_color_map,
        label_kws=dict(size=14, color="black", r=110, orientation="vertical"),
        link_kws_handler=link_handler,
    )
    return circos


def plot_circos_plot(df: pd.DataFrame, color="aggregated_interaction", color_map: Dict[str, str] | None = None, niche_name: str | None = None, aggregator: str | None = None):
    assert color in ["aggregated_interaction", "above_median_fraction"]
    assert "above_median_fraction" and "aggregated_interaction" in df.columns
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    circos = get_circos(df=df, color=color, color_map=color_map)
    fig = circos.plotfig(figsize=(20, 15))
    fig.subplots_adjust(top=3, right=3)

    niche_name = niche_name or ""
    aggregator = f"({aggregator}) " if aggregator else ""
    fig.suptitle(
        f"Circos Plot of Aggregated Interactions {aggregator}in Niche {niche_name}\n"
        f"Color = {color}, Line Width = {[col for col in df.columns.tolist() if col != color][0]}",
        fontsize=16,
        fontweight="bold",
        y=0.99,
    )

    val_min = min(df[color])
    val_max = max(df[color])
    norm = colors.Normalize(vmin=val_min, vmax=val_max)
    cmap = cm.get_cmap("Reds")
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar_ax = fig.add_axes([0.91, 0.25, 0.02, 0.5])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical", label=color)
    cbar.set_label(color, fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    return fig


def main(data_dir: Path | None = None):
    data_dir = data_dir or resolve_data_dir()
    save_dir = resolve_output_figures_dir() / "figure7"
    save_dir.mkdir(parents=True, exist_ok=True)
    resources_dir = Path(__file__).resolve().parents[2] / "resources"

    with open(resources_dir / "colormaps.yaml") as f:
        colormaps = yaml.safe_load(f)
    color_dict_label = colormaps["label"]

    interactions_dir = data_dir / "niches" / "interactions"

    for panel, niche in NICHES.items():
        df_path = interactions_dir / f"{niche}.parquet"
        df_example = pd.read_parquet(df_path)
        df_example = df_example.dropna()

        df_circos = df_example.set_index(["Cell_Type_1", "Cell_Type_2"])
        df_circos.rename(columns={"Interaction_Frequency": "aggregated_interaction", "Above_Median_Fraction": "above_median_fraction"}, inplace=True)
        df_circos = df_circos[["aggregated_interaction", "above_median_fraction"]]

        fig = plot_circos_plot(df=df_circos, color="above_median_fraction", color_map=color_dict_label, niche_name=niche)
        fig_path = save_dir / f"figure7{panel}_{niche}_circos_plot.pdf"
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        logger.info(f"saved Figure 7{panel} (niche {niche}) to {fig_path}")


if __name__ == "__main__":
    CLI(main)
