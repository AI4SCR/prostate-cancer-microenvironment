"""Reproduce Supplementary Figure 5a: stacked barplot of niche composition
per TMA core, for cores with high `tumorERG+p53+_ProlifLuminal` (niche 6)
abundance, with pat_id/disease_progr/gleason_grp/inflammation annotation
rows on top.

1:1 port of `sync_paper/06-spatial-niches/abundance/stacked_frequencies.py`
(the `group_var == 'tma_id'` branch, lines ~89-251 of that file) --
identified by Melissa (paper co-author) as the actual source for this
panel; a prior pass over this file wrongly ruled it out based on a
mis-transcribed legend (see `figures.md`'s S5a Issue note).

Disclosed extension, per Melissa's direct guidance ("Andrea added some
metadata annotation manually... I think it should be quite easy to also
code this"): the source's `ann_rows = ["disease_progr", "gleason_grp"]`
is genuinely broken as literally written -- `tma_cols` (the column subset
kept for the `tma_id` branch) never includes `disease_progr`, so
`df_metadata_p53["disease_progr"]` would `KeyError`. Extended `tma_cols`
to include `disease_progr`, and added a `pat_id` column (merged in from
`clinical.parquet` via `tma_id`, since `pat_id` is only kept on the
`pat_id`-branch's column list in the source) and `pat_id`/`inflammation`
to `ann_rows`, matching the 4-row annotation (PatID, Disease progr,
Gleason grp, Inflammation) visible in the published panel. Colors for both
come from `resources/colormaps.yaml`'s existing `pat_id`/`disease_progr`
sections (already used elsewhere in this repo) -- not invented.

Second disclosed deviation, for the rebuttal letter's Issue 2 response:
`THRESHOLD` was raised from the verbatim source's `0.05` to `0.1195737`
so this panel selects the same 5 patients as Fig 6b's niche-6 KM plot
(`figure6_km_niche6.R`), rather than 6. The two panels used different
"positive for niche 6" cutoffs -- S5a's per-core frequency (with
`pseudocount=1`) at 0.05 vs. Fig 6b's median split over nonzero cores'
cell-level frequency (`pseudocount=0`) -- which meant S5a included one
extra, non-progressing patient (`pat_id 96.22128`) that Fig 6b's stricter
cutoff excludes. `0.1195737` is Fig 6b's own median threshold value,
verified (via an ad hoc check, not part of this script) to select the
same 8 TMA cores / 5 patients under either frequency formula at this
cutoff. Confirmed empirically: at the new threshold, all 5 selected
patients progressed (`95.20582, 96.7481, 97.1247, 97.5521, 98.6114`).
"""

from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from jsonargparse import CLI
from loguru import logger

from prostate_cancer.utils import resolve_data_dir, resolve_output_figures_dir

NICHE_OF_INTEREST = "tumorERG+p53+_ProlifLuminal"  # niche 6
THRESHOLD = 0.1195737  # matches Fig 6b's median-split cutoff; was 0.05 -- see module docstring


def compute_label_frequency(data: pd.DataFrame, level: str, pseudocount: int = 1, group_vars: list[str] = ["sample_id"]) -> pd.Series:
    data[level] = data[level].astype("category")
    if pseudocount > 0:
        pdat = data.groupby(group_vars, observed=False)[level].value_counts()
        pdat += 1
        pdat /= pdat.groupby(group_vars, observed=False).sum()
        pdat.name = "proportion"
    else:
        pdat = data.groupby(group_vars, observed=False)[level].value_counts(normalize=True)
    return pdat


def get_label_frequency_table(data: pd.DataFrame, level: str, group_vars: list[str] = ["sample_id"]) -> pd.DataFrame:
    props = compute_label_frequency(data=data, level=level, pseudocount=1, group_vars=group_vars)
    props = props.reset_index().pivot(index=group_vars, columns=level, values="proportion")
    props.columns = props.columns.astype(str)
    return props.astype(float)


def main(data_dir: Path | None = None):
    data_dir = data_dir or resolve_data_dir()
    figures_dir = resolve_output_figures_dir() / "figureS5"
    figures_dir.mkdir(parents=True, exist_ok=True)
    resources_dir = Path(__file__).resolve().parents[2] / "resources"

    var_name = "niche"
    group_var = "tma_id"

    df_clusters = pd.read_parquet(data_dir / "niches" / "clusters_annotated.parquet")
    logger.info(f"clusters shape: {df_clusters.shape}")

    clinical = pd.read_parquet(data_dir / "clinical.parquet")

    valid_tma_ids = df_clusters["tma_id"].unique().tolist()
    df_metadata = clinical[clinical["tma_id"].isin(valid_tma_ids)].copy()

    # pat_id map, per Melissa's guidance (not in tma_cols on the source's tma_id branch)
    pat_id_map = df_metadata.drop_duplicates(subset=["tma_id"]).set_index("tma_id")["pat_id"]

    df_metadata.set_index("tma_id", inplace=True)

    # tma_cols: source list + disclosed addition of "disease_progr" (source's own
    # ann_rows references it, but it's missing from tma_cols as literally written)
    tma_cols = ["stromogenic_smc_loss_reactive_stroma_present", "inflammation", "glandular_atrophy_pin", "gleason_grp", "disease_progr"]

    metadata = df_metadata[tma_cols].copy()
    metadata.reset_index(inplace=True)
    metadata.drop_duplicates(subset=[group_var], inplace=True)
    metadata.set_index(group_var, inplace=True)
    metadata["pat_id"] = pat_id_map.reindex(metadata.index)

    df_freqs = get_label_frequency_table(data=df_clusters, level=var_name, group_vars=[group_var])
    df_freqs, metadata = df_freqs.align(metadata, join="inner", axis=0)
    logger.info(f"frequencies shape after alignment: {df_freqs.shape}")
    logger.info(f"metadata shape after alignment: {metadata.shape}")

    # make long
    df_freqs_long = df_freqs.reset_index().melt(id_vars=group_var, var_name=var_name, value_name="frequency")
    df_freqs_long = df_freqs_long.merge(metadata.reset_index(), on=group_var, how="left")

    niche_of_interest = NICHE_OF_INTEREST
    threshold = THRESHOLD

    # samples to include (high for niche_of_interest)
    df_freqs_p53 = df_freqs_long.loc[df_freqs_long[var_name] == niche_of_interest].copy()
    df_freqs_p53["high_freq"] = df_freqs_p53["frequency"] > threshold
    df_freqs_p53_high = df_freqs_p53.loc[df_freqs_p53["high_freq"]].copy()

    high_samples = df_freqs_p53_high[group_var].unique()
    df_long_p53 = df_freqs_long.loc[df_freqs_long[group_var].isin(high_samples)].copy()

    # order TMAs by descending frequency of niche_of_interest
    tma_order = df_freqs_p53_high.sort_values("frequency", ascending=False)[group_var].tolist()

    # niches that should come first in the stack
    niches_first = [
        "tumorERG+p53+_ProlifLuminal",
        "Macrophages_Tcells_CAF1(CD105-)",
        "immune_bloodvessels_CAF1(CD105-)",
        "TLS",
    ]

    # IMPORTANT: aggregate in case you have repeated rows per (sample, niche)
    df_plot = df_long_p53.groupby([group_var, var_name], as_index=False)["frequency"].sum()

    # build full niche order: niches_first, then all remaining niches (sorted by overall abundance)
    all_niches = df_plot[var_name].unique().tolist()
    rest = [n for n in all_niches if n not in niches_first]

    rest_totals = df_plot.loc[df_plot[var_name].isin(rest)].groupby(var_name)["frequency"].sum().sort_values(ascending=False)
    rest = rest_totals.index.tolist()

    niche_order = [n for n in niches_first if n in all_niches] + rest

    # colors
    with open(resources_dir / "colormaps.yaml") as f:
        colormaps = yaml.safe_load(f)
    color_dict_niche = colormaps["niche"]

    # pivot wide for stacking
    df_wide = df_plot.pivot(index=group_var, columns=var_name, values="frequency").reindex(index=tma_order, columns=niche_order).fillna(0.0)

    # --- prepare metadata in correct order
    df_metadata_p53 = metadata.loc[high_samples].copy()
    df_metadata_p53 = df_metadata_p53.reindex(index=tma_order)

    # --- figure with annotation panel + stacked barplot
    fig = plt.figure(figsize=(14, 7))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[0.6, 6], hspace=0.05)

    ax_ann = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharex=ax_ann)

    # -------------------------
    # stacked bar plot
    # -------------------------
    bottom = np.zeros(df_wide.shape[0])
    x = np.arange(df_wide.shape[0])

    for niche in niche_order:
        heights = df_wide[niche].values
        ax.bar(
            x,
            heights,
            bottom=bottom,
            label=niche,
            color=color_dict_niche.get(niche, "#333333"),
            width=1.0,  # bars touch
            edgecolor="black",  # add border
            linewidth=0.4,  # thin border looks better
        )
        bottom += heights

    # labels
    ax.set_xlabel("TMA ID")
    ax.set_ylabel("Frequency")

    # remove title
    ax.set_title(None)

    # ticks
    ax.set_xticks(x)
    ax.set_xticklabels(df_wide.index, rotation=45, ha="right")

    # remove padding so bars touch edges
    ax.set_xlim(-0.5, len(x) - 0.5)

    # remove box
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

    # legend
    ax.legend(title="Niche", bbox_to_anchor=(1.02, 1), loc="upper left", ncol=1)

    # -------------------------
    # annotation panel on top
    # -------------------------
    ann_rows = ["pat_id", "disease_progr", "gleason_grp", "inflammation"]
    ann_color_dicts = {
        "pat_id": colormaps["pat_id"],
        "disease_progr": colormaps["disease_progr"],
        "gleason_grp": colormaps["gleason_grp"],
        # yaml.safe_load() parses bareword yes/no keys as YAML 1.1 booleans
        # (True/False), not the strings "yes"/"no" the data actually uses --
        # same gotcha found and fixed in figure6_niche_abundance_heatmap.R
        "inflammation": {("yes" if k is True else "no" if k is False else k): v for k, v in colormaps["inflammation"].items()},
    }

    bar_width = 1.0

    for row_idx, ann_var in enumerate(ann_rows):
        y = len(ann_rows) - 1 - row_idx  # puts first variable on top
        vals = df_metadata_p53[ann_var]

        for i, v in enumerate(vals):
            # convert to string for dictionary lookup
            key = str(v)

            color = ann_color_dicts[ann_var].get(key, "#d3d3d3")

            rect = patches.Rectangle(
                (i - bar_width / 2, y),
                bar_width,
                1,
                facecolor=color,
                edgecolor="white",
                linewidth=0.5,
            )
            ax_ann.add_patch(rect)

    # annotation axis styling
    ax_ann.set_xlim(-0.5, len(x) - 0.5)
    ax_ann.set_ylim(0, len(ann_rows))
    ax_ann.set_yticks(np.arange(len(ann_rows)) + 0.5)
    ax_ann.set_yticklabels(ann_rows[::-1])
    ax_ann.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ["top", "right", "bottom", "left"]:
        ax_ann.spines[spine].set_visible(False)

    fig.tight_layout()  # adjust layout to prevent overlap
    fig_path = figures_dir / "p53_immune_full_barplot_annotated.pdf"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    logger.info(f"saved Supplementary Figure 5a to {fig_path}")


if __name__ == "__main__":
    CLI(main)
