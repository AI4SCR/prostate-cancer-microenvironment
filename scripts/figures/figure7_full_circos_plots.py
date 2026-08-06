"""1:1 port of `sync_paper/06-spatial-niches/interactions/interaction_compute_circos.py`.

Unlike `figure7def_circos_plots.py` (ported from the older
`000_paper/11_niches/114_interactions/visualize_circos_plot.py` +
`circos_plots.py`, which only plots precomputed `dataframes_v2/*.parquet`),
this is the newer `sync_paper` source's *full* pipeline: it computes the
per-sample, per-niche cell-cell interaction matrices itself from the raw
anndata pickles and `clusters_annotated.parquet` (not `_v2`), then plots the
same three niches. Both scripts are kept -- see `figures.md`.

Only path/env-var substitution and one disclosed environment-compatibility
patch (see `_patch_anndata_unpickling` below) are applied; no computational
logic is changed.
"""

import pickle
import weakref
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
from tqdm import tqdm

from prostate_cancer.utils import resolve_legacy_dir, resolve_output_figures_dir

matplotlib.use("Agg")

# Live legacy path (not staged under LEGACY_DATA_DIR -- 476 pickles, too
# large to duplicate into the repo); matches legacy's own hardcoded
# `.../PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas`
# path exactly, confirmed 100% sample_id overlap with clusters_annotated_v2.parquet
# (see data/assets.md). Read live per CLAUDE.md's "(b) legacy data ... read
# live from a legacy path" allowance.
ANNDATAS_DIR = Path(
    "/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment"
    "/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas"
)

NICHES_WITH_INTERACTION = ["luminal_infiltrated", "luminal_CAF1(CD105High)", "tumor_CAF1(CD105High)"]


def _patch_anndata_unpickling():
    """Environment-compatibility patch, not a logic change.

    These pickles were written by an older `anndata` whose `AnnDataFileManager`
    stored its backing-file reference under the state key `_adata`; the
    installed `anndata` (0.12.10) expects `_adata_ref` and raises `KeyError`
    otherwise. Confirmed these objects aren't actually file-backed
    (`_filename`/`_file` are `None`), so this is a pure key-rename shim.
    """
    from anndata._core.file_backing import AnnDataFileManager

    def patched_setstate(self, state):
        if "_adata_ref" not in state and "_adata" in state:
            state["_adata_ref"] = state.pop("_adata")
        self.__dict__["_adata_ref"] = weakref.ref(state["_adata_ref"])
        self.__dict__["_filename"] = state.get("_filename")
        self.__dict__["_filemode"] = state.get("_filemode")
        self.__dict__["_file"] = state.get("_file")

    AnnDataFileManager.__setstate__ = patched_setstate


def compute_interaction_matrices(adata, cell_types=None, graph_key="radius_32"):
    assert adata.n_obs > 15, "Anndata object has too few observations to compute interactions."
    if cell_types is None:
        cell_types = adata.obs["label"].unique().tolist()
        cell_types.sort()
    interaction_matrix = pd.DataFrame(0, index=cell_types, columns=cell_types)
    adj = adata.obsp[graph_key].toarray()
    np.fill_diagonal(adj, 0)
    from scipy.sparse import csr_matrix

    sparse_adj = csr_matrix(adj)
    adata.obsp[f"{graph_key}_no_self_loops"] = sparse_adj
    from athena.utils.general import get_nx_graph_from_anndata

    g = get_nx_graph_from_anndata(adata, key=f"{graph_key}_no_self_loops")
    edges = list(g.edges())
    sample_interaction = interaction_matrix.copy()
    for edge in edges:
        (c1, c2) = edge
        (type1, type2) = (adata.obs.loc[c1]["label"], adata.obs.loc[c2]["label"])
        if type1 > type2:
            (type1, type2) = (type2, type1)
        sample_interaction.loc[type1, type2] += 1
    num_edges = len(edges)
    print(f"Total number of edges in graph: {num_edges}")
    sample_interaction_freqs = sample_interaction / num_edges
    sample_interaction_freqs = sample_interaction_freqs.T
    mask = np.triu(np.ones_like(sample_interaction_freqs, dtype=bool), k=1)
    sample_interaction_freqs = sample_interaction_freqs.mask(mask)
    return sample_interaction_freqs


def plot_interaction_matrix(interaction_freqs, title=None):
    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(16, 12))
    sns.heatmap(interaction_freqs, annot=True, fmt=".2f", cmap="Reds", cbar_kws={"label": "Interaction Frequency"})
    plt.title(title)
    plt.xlabel("Cell Type 2")
    plt.ylabel("Cell Type 1")
    plt.tight_layout()
    return plt.gcf()


def check_triangular_zero_xor(df: pd.DataFrame):
    if df.shape[0] != df.shape[1]:
        return False, "DataFrame is not square. Cannot perform triangular check."
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
    assert color in ["aggregated_interaction", "above_median_fraction"], "color must be either 'aggregated_interaction' or 'above_median_fraction'"
    assert "above_median_fraction" and "aggregated_interaction" in df.columns
    assert df.index.nlevels == 2, "DataFrame must have a MultiIndex "

    if df.index.names != ["cell_type_1", "cell_type_2"]:
        df.index.set_names(["cell_type_1", "cell_type_2"], inplace=True)

    width = [col for col in df.columns.tolist() if col != color][0]
    width_df = df[width]
    width_df = width_df.reset_index()
    width_df.columns = ["from", "to", "Value"]
    width_df["Value"].dropna()
    width_dict = {(row["from"], row["to"]): row["Value"] for _, row in width_df.iterrows()}

    color_df = df[color]
    color_df = color_df.reset_index()
    color_df.columns = ["from", "to", "Value"]
    color_df["Value"].dropna()
    color_dict = {(row["from"], row["to"]): row["Value"] for _, row in color_df.iterrows()}

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


def plot_circos_plot(df: pd.DataFrame, color="aggregated_interaction", color_map=None, niche_name=None, aggregator=None, save_path=None):
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

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def main(legacy_dir: Path | None = None):
    _patch_anndata_unpickling()

    legacy_dir = legacy_dir or resolve_legacy_dir()
    output_dir = resolve_output_figures_dir() / "figure7" / "full_circos" / "data"
    figures_dir = resolve_output_figures_dir() / "figure7" / "full_circos"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    resources_dir = Path(__file__).resolve().parents[2] / "resources"
    logger.info(f"Results will be saved to {output_dir}")
    logger.info(f"Figures will be saved to {figures_dir}")

    ####### CELLTYPE ENRICHMENT HEATMAPS #######
    df_clusters = pd.read_parquet(legacy_dir / "5-niches" / "annotation" / "clusters_annotated.parquet")
    logger.info(f"Clusters shape: {df_clusters.shape}")
    # legacy bug: reset_index(inplace=True) here destroyed the sample_id/
    # object_id MultiIndex that the very next line (get_level_values) and
    # the later df_clusters.loc[sample] both require -- removed, retaining
    # the index clusters_annotated.parquet is already stored with; no other
    # logic changed (see bugs.md)

    anndata_dir = ANNDATAS_DIR

    sample_list = df_clusters.index.get_level_values("sample_id").unique().tolist()

    ############ compute interactions ##################
    cell_types = df_clusters["label"].unique().tolist()
    cell_types.sort()
    graph_key = "radius_32"

    niche_ids = df_clusters["niche"].unique().tolist()
    niche_ids.sort()

    anndata_save_dir = output_dir / "anndatas_with_interactions"
    anndata_save_dir.mkdir(parents=True, exist_ok=True)

    for sample in tqdm(sample_list):
        logger.info(f"Processing sample: {sample}")
        anndata_path = anndata_dir / f"{sample}.pkl"
        if not anndata_path.exists():
            logger.info(f"Anndata file {anndata_path} does not exist")
        with open(anndata_path, "rb") as f:
            adata = pickle.load(f)
        obs = adata.obs
        sample_data = df_clusters.loc[sample]
        logger.info(f"shape of sample data: {sample_data.shape} ")
        logger.info(f"shape of obs: {obs.shape} ")
        sample_data.index = sample_data.index.astype(str)

        obs = obs.merge(sample_data, left_index=True, right_index=True, how="left", suffixes=("", "_full"))
        adata.obs = obs

        interaction_freqs = compute_interaction_matrices(adata, cell_types=cell_types, graph_key=graph_key)
        title = f"Cell-Cell Interaction Matrix for Sample {sample}"
        plot_interaction_matrix(interaction_freqs, title=title)

        adata.uns["interaction_freqs_overall"] = interaction_freqs

        for niche in niche_ids:
            logger.info(f"Processing niche: {niche}")
            niche_adata = adata[adata.obs["niche"] == niche].copy()
            if niche_adata.n_obs <= 15:
                logger.info(f"Skipping niche {niche} with only {niche_adata.n_obs} cells")
                continue
            interaction_freqs = compute_interaction_matrices(niche_adata, cell_types=cell_types, graph_key=graph_key)
            adata.uns[f"interaction_freqs_niche_{niche}"] = interaction_freqs
            title = f"Cell-Cell Interaction Matrix for Sample {sample}, Niche {niche}"
            plot_interaction_matrix(interaction_freqs, title=title)

        anndata_save_path = anndata_save_dir / f"{sample}_with_interactions.pkl"
        anndata_save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(anndata_save_path, "wb") as f:
            pickle.dump(adata, f)

    ############ get interaction matrices ##################
    cell_types = df_clusters["label"].unique().tolist()
    cell_types.sort()
    graph_key = "radius_32"

    niche_ids = df_clusters["niche"].unique().tolist()
    niche_ids.sort()

    anndata_dir = anndata_save_dir
    overall_interactions = []
    overall_key = "interaction_freqs_overall"

    for sample_id in tqdm(sample_list):
        logger.info(f"  Sample: {sample_id}")
        anndata_path = anndata_dir / f"{sample_id}_with_interactions.pkl"

        with open(anndata_path, "rb") as f:
            adata = pickle.load(f)

        if overall_key in adata.uns:
            interaction_freqs = adata.uns[overall_key]
            interaction_df_long = interaction_freqs.reset_index().melt(id_vars="index")
            interaction_df_long.columns = ["Cell_Type_1", "Cell_Type_2", "Interaction_Frequency"]
            interaction_df_long["Sample_ID"] = sample_id
            overall_interactions.append(interaction_df_long)
        else:
            logger.info(f"Interaction frequencies not found in sample {sample_id}")

    overall_interaction_freqs = pd.concat(overall_interactions, ignore_index=True)
    overall_interaction_freqs["Interaction_Frequency"] = overall_interaction_freqs["Interaction_Frequency"].fillna(0)
    overall_medians = overall_interaction_freqs.groupby(["Cell_Type_1", "Cell_Type_2"])["Interaction_Frequency"].median()
    overall_interaction_freqs.set_index(["Cell_Type_1", "Cell_Type_2"], inplace=True)

    dict_results = {}
    for niche_id in tqdm(niche_ids):
        logger.info(f"Processing niche: {niche_id}")
        key_name = f"interaction_freqs_niche_{niche_id}"

        interactions = []
        above_median = pd.DataFrame(index=overall_medians.index, columns=["Above_Median"])
        above_median["Above_Median"] = 0
        below_median = pd.DataFrame(index=overall_medians.index, columns=["Below_Median"])
        below_median["Below_Median"] = 0
        num_samples_with_niche = 0
        for sample_id in sample_list:
            logger.info(f"  Sample: {sample_id}")
            anndata_path = anndata_dir / f"{sample_id}_with_interactions.pkl"

            with open(anndata_path, "rb") as f:
                adata = pickle.load(f)

            if key_name in adata.uns:
                interaction_freqs = adata.uns[key_name]
                interaction_df_long = interaction_freqs.reset_index().melt(id_vars="index")
                interaction_df_long.columns = ["Cell_Type_1", "Cell_Type_2", "Interaction_Frequency"]
                interaction_df_long.set_index(["Cell_Type_1", "Cell_Type_2"], inplace=True)
                assert overall_medians.shape[0] == interaction_df_long.shape[0], "Mismatch in interaction frequencies shape."
                medians_aligned, interaction_df_long = overall_medians.align(interaction_df_long, join="inner", axis=0)
                interaction_df_long["lfc"] = np.log2((interaction_df_long["Interaction_Frequency"] + 1e-9) / (medians_aligned + 1e-9))
                comparison = interaction_df_long["Interaction_Frequency"] > medians_aligned
                comparison_negative = interaction_df_long["Interaction_Frequency"] < medians_aligned
                above_median.loc[comparison.index, "Above_Median"] += comparison.astype(int)
                below_median.loc[comparison_negative.index, "Below_Median"] += comparison_negative.astype(int)
                interactions.append(interaction_df_long)
                num_samples_with_niche += 1
            else:
                logger.info(f"Interaction frequencies for niche {niche_id} not found in sample {sample_id}")

        above_median["Above_Median_Fraction"] = above_median["Above_Median"] / num_samples_with_niche
        below_median["Below_Median_Fraction"] = below_median["Below_Median"] / num_samples_with_niche

        mean_interactions = sum(interactions) / len(interactions)

        df_plot = mean_interactions.merge(above_median, left_index=True, right_index=True).merge(below_median, left_index=True, right_index=True)
        df_plot = df_plot.reset_index()
        df_plot["size"] = np.where(df_plot["lfc"] > 0, df_plot["Above_Median_Fraction"], df_plot["Below_Median_Fraction"])
        dict_results[niche_id] = df_plot

    results_dir = output_dir / "dataframes"
    results_dir.mkdir(parents=True, exist_ok=True)
    for niche_id, df_result in dict_results.items():
        df_path = results_dir / f"{niche_id}.parquet"
        df_result.to_parquet(df_path)

    ### plot circos plots for each niche of interest
    data_dir = results_dir

    with open(resources_dir / "colormaps.yaml") as f:
        color_maps = yaml.safe_load(f)
    color_dict_label = color_maps["label"]

    plot_dir = figures_dir / "circos_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for niche in NICHES_WITH_INTERACTION:
        df_path = data_dir / f"{niche}.parquet"
        df_example = pd.read_parquet(df_path)
        df_example = df_example.dropna()

        df_circos = df_example.set_index(["Cell_Type_1", "Cell_Type_2"])
        df_circos.rename(columns={"Interaction_Frequency": "aggregated_interaction", "Above_Median_Fraction": "above_median_fraction"}, inplace=True)
        df_circos = df_circos[["aggregated_interaction", "above_median_fraction"]]
        fig_path = plot_dir / f"{niche}_circos_plot.pdf"
        plot_circos_plot(df=df_circos, color="above_median_fraction", color_map=color_dict_label, niche_name=niche, aggregator=None, save_path=fig_path)
        logger.info(f"saved circos plot for niche {niche} to {fig_path}")


if __name__ == "__main__":
    CLI(main)
