"""
RQ3: Immune Neighborhood Composition and TLS-like Structure Characterization

Builds radius-32 spatial graphs per ROI, profiles the local cellular neighborhoods
of B cells, computes a TLS-like co-localization score, and tests associations with
histological inflammation annotation and PSA/clinical outcome.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ3")
OUT_TAB = Path("output/tables/RQ3")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load data ─────────────────────────────────────────────────────────────

print("Loading data...")
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
spatial_files = sorted((BASE_DIR / "02_processed/features/spatial/filtered-annotated").glob("*.parquet"))

# Build lookup by sample_id stem
meta_by_id = {f.stem: f for f in meta_files}
spatial_by_id = {f.stem: f for f in spatial_files}

common_ids = sorted(set(meta_by_id) & set(spatial_by_id))
assert len(common_ids) > 0

clinical = pd.read_parquet(BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet")
tumor_ids = set(clinical.index[clinical.is_tumor == "yes"])

# Filter to tumor ROIs
sample_ids = [sid for sid in common_ids if sid in tumor_ids]
print(f"  Processing {len(sample_ids)} tumor ROIs")

# %% ─── Cell type definitions ─────────────────────────────────────────────────

B_CELL_LABEL = "immune-B-cells(CD20+)"
T_HELPER_LABEL = "immune-T-cells_helper(CD3+CD4+)"
T_REG_LABEL = "immune-T-cells_regulatory(CD3+CD4+FoxP3+)"
T_CYTOTOXIC_LABEL = "immune-T-cells_cytotoxic(CD3+CD8a+)"
T_CELL_LABELS = {T_HELPER_LABEL, T_REG_LABEL, T_CYTOTOXIC_LABEL, "immune-T-cells(CD3+)"}
IMMUNE_LABELS = {B_CELL_LABEL, T_HELPER_LABEL, T_REG_LABEL, T_CYTOTOXIC_LABEL,
                 "immune-T-cells(CD3+)", "immune-macrophages(CD68+)",
                 "immune-T-helper-B-cells", "immune-T-helper-T-cytotoxic",
                 "immune-T-helper-macrophages", "immune-PMN-MDSCs(CD11b+CD66b+)"}

RADIUS = 32  # pixels, as in draft paper

# %% ─── Per-ROI neighborhood analysis ────────────────────────────────────────

print("Computing B-cell neighborhood profiles per ROI...")


def compute_bcell_neighborhood_profile(spatial_df: pd.DataFrame, meta_df: pd.DataFrame, radius: float) -> dict:
    """
    For each B cell in an ROI, extract the label composition of its radius neighborhood.
    Returns aggregate statistics for the ROI.
    """
    # Align indices
    common_obj = spatial_df.index.intersection(meta_df.index)
    if len(common_obj) == 0:
        return {}

    spatial_df = spatial_df.loc[common_obj]
    meta_df = meta_df.loc[common_obj]

    coords = spatial_df[["x", "y"]].values
    labels = meta_df["label"].values

    b_mask = labels == B_CELL_LABEL
    n_bcells = b_mask.sum()
    n_cells = len(labels)

    if n_bcells == 0:
        return {"n_bcells": 0, "n_cells": n_cells}

    # Build KD-tree
    tree = cKDTree(coords)

    # For each B cell, find neighbors within radius
    b_coords = coords[b_mask]
    neighbor_lists = tree.query_ball_point(b_coords, r=radius)

    # Count label frequencies in neighborhoods (excluding the B cell itself)
    b_indices = np.where(b_mask)[0]
    neighbor_label_counts = {}
    t_cell_coloc_count = 0  # B cells with at least one T cell neighbor
    t_helper_coloc_count = 0
    t_reg_coloc_count = 0

    for i, (b_idx, neighbors) in enumerate(zip(b_indices, neighbor_lists)):
        neighbor_labels = labels[neighbors]
        neighbor_labels = neighbor_labels[np.array(neighbors) != b_idx]  # exclude self

        has_tcell = any(l in T_CELL_LABELS for l in neighbor_labels)
        if has_tcell:
            t_cell_coloc_count += 1
        has_thelper = any(l == T_HELPER_LABEL for l in neighbor_labels)
        if has_thelper:
            t_helper_coloc_count += 1
        has_treg = any(l == T_REG_LABEL for l in neighbor_labels)
        if has_treg:
            t_reg_coloc_count += 1

        for lbl in neighbor_labels:
            neighbor_label_counts[lbl] = neighbor_label_counts.get(lbl, 0) + 1

    # Normalize by total neighbor count
    total_neighbor_hits = sum(neighbor_label_counts.values())

    result = {
        "n_bcells": n_bcells,
        "n_cells": n_cells,
        "bcell_proportion": n_bcells / n_cells,
        "tls_score": t_cell_coloc_count / n_bcells,  # fraction of B cells with T cell neighbor
        "tls_thelper_score": t_helper_coloc_count / n_bcells,
        "tls_treg_score": t_reg_coloc_count / n_bcells,
    }

    # Add neighborhood label fractions
    for lbl, count in neighbor_label_counts.items():
        safe_lbl = lbl.replace("(", "_").replace(")", "").replace("+", "pos").replace("-", "neg").replace(" ", "_")
        result[f"bcell_nbhd_{safe_lbl}"] = count / max(total_neighbor_hits, 1)

    return result


roi_profiles = []
for sid in sample_ids:
    spatial_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")

    # Strip sample_id from index for easy alignment
    spatial_df = spatial_df.loc[sid] if sid in spatial_df.index.get_level_values("sample_id") else spatial_df
    meta_df = meta_df.loc[sid] if sid in meta_df.index.get_level_values("sample_id") else meta_df

    profile = compute_bcell_neighborhood_profile(spatial_df, meta_df, radius=RADIUS)
    profile["sample_id"] = sid
    roi_profiles.append(profile)

profiles_df = pd.DataFrame(roi_profiles).set_index("sample_id")
profiles_df = profiles_df.fillna(0)
print(f"  Computed profiles for {len(profiles_df)} ROIs")
print(f"  ROIs with ≥1 B cell: {(profiles_df.n_bcells > 0).sum()}")
print(f"  Mean TLS score (B cells with T cell neighbor): {profiles_df['tls_score'].mean():.3f}")

profiles_df.to_csv(OUT_TAB / "roi_bcell_neighborhood_profiles.csv")

# %% ─── Immunotype annotation ─────────────────────────────────────────────────

# Merge with clinical annotations
profiles_with_clin = profiles_df.join(clinical[["inflammation", "psa_progr", "psa_progr_time",
                                                  "clinical_progr", "clinical_progr_time",
                                                  "gleason_grp", "pat_id"]], how="left")

# Only ROIs with ≥3 B cells for TLS analysis (to reduce noise)
bcell_rich = profiles_with_clin[profiles_with_clin.n_bcells >= 3].copy()
print(f"\nROIs with ≥3 B cells: {len(bcell_rich)}")

# %% ─── Figure 1: TLS score vs inflammation ───────────────────────────────────

print("Generating figures...")

# 1a. TLS score by inflammation status (ROIs with ≥1 B cell)
bcell_present = profiles_with_clin[profiles_with_clin.n_bcells >= 1].dropna(subset=["inflammation"])
infl_yes = bcell_present[bcell_present.inflammation == "yes"]["tls_score"].values
infl_no = bcell_present[bcell_present.inflammation == "no"]["tls_score"].values

stat, pval = mannwhitneyu(infl_yes, infl_no, alternative="two-sided")
print(f"  TLS score inflamed vs non-inflamed: p={pval:.4f}, n_yes={len(infl_yes)}, n_no={len(infl_no)}")

fig, axes = plt.subplots(1, 3, figsize=(14, 5))

ax = axes[0]
sns.boxplot(data=bcell_present, x="inflammation", y="tls_score", order=["no", "yes"],
            ax=ax, palette={"no": "#aec6cf", "yes": "#ff7f7f"})
ax.set_title(f"TLS score by inflammation\np={pval:.3e}")
ax.set_xlabel("Histological inflammation")
ax.set_ylabel("TLS score (fraction of B cells\nwith ≥1 T cell neighbor)")

# 1b. B-cell proportion by inflammation
bcell_all = profiles_with_clin.dropna(subset=["inflammation"])
stat2, pval2 = mannwhitneyu(
    bcell_all[bcell_all.inflammation == "yes"]["bcell_proportion"].values,
    bcell_all[bcell_all.inflammation == "no"]["bcell_proportion"].values,
    alternative="two-sided"
)
ax = axes[1]
sns.boxplot(data=bcell_all, x="inflammation", y="bcell_proportion", order=["no", "yes"],
            ax=ax, palette={"no": "#aec6cf", "yes": "#ff7f7f"})
ax.set_title(f"B-cell proportion by inflammation\np={pval2:.3e}")
ax.set_xlabel("Histological inflammation")
ax.set_ylabel("B-cell proportion")

# 1c. TLS Treg score vs T helper score by inflammation
ax = axes[2]
colors = {"yes": "#ff7f7f", "no": "#aec6cf"}
for inf_val in ["yes", "no"]:
    sub = bcell_present[bcell_present.inflammation == inf_val]
    ax.scatter(sub["tls_thelper_score"], sub["tls_treg_score"],
               c=colors[inf_val], alpha=0.5, s=20, label=f"inflammation={inf_val}")
ax.set_xlabel("T helper coloc score")
ax.set_ylabel("T regulatory coloc score")
ax.set_title("T helper vs Treg co-localization\nwith B cells")
ax.legend()

plt.tight_layout()
fig.savefig(OUT_FIG / "TLS_score_inflammation.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 2: Neighborhood composition heatmap ────────────────────────────

# Get neighborhood label columns
nbhd_cols = [c for c in profiles_df.columns if c.startswith("bcell_nbhd_")]

# Mean neighborhood composition: inflamed vs non-inflamed
bcell_nbhd = profiles_with_clin[nbhd_cols + ["inflammation"]].dropna(subset=["inflammation"])
mean_nbhd = bcell_nbhd.groupby("inflammation")[nbhd_cols].mean()

# Rename columns for display
col_rename = {c: c.replace("bcell_nbhd_", "").replace("_", " ") for c in nbhd_cols}
mean_nbhd = mean_nbhd.rename(columns=col_rename)

# Filter to top columns
top_cols = mean_nbhd.max().sort_values(ascending=False).head(15).index.tolist()
mean_nbhd_top = mean_nbhd[top_cols]

fig, ax = plt.subplots(figsize=(12, 4))
sns.heatmap(
    mean_nbhd_top.T,
    cmap="YlOrRd",
    ax=ax,
    cbar_kws={"label": "Mean fraction in B-cell neighborhood"},
    annot=True,
    fmt=".3f",
    annot_kws={"size": 8},
)
ax.set_title("B-cell neighborhood composition by inflammation status")
ax.set_xlabel("Inflammation status")
ax.set_ylabel("Cell type in neighborhood")
plt.tight_layout()
fig.savefig(OUT_FIG / "bcell_neighborhood_heatmap.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 3: TLS score vs PSA recurrence ─────────────────────────────────

bcell_psa = profiles_with_clin[profiles_with_clin.n_bcells >= 1].dropna(subset=["psa_progr"])
stat3, pval3 = mannwhitneyu(
    bcell_psa[bcell_psa.psa_progr == 1]["tls_score"].values,
    bcell_psa[bcell_psa.psa_progr == 0]["tls_score"].values,
    alternative="two-sided"
)
print(f"  TLS score vs PSA recurrence: p={pval3:.4f}")

bcell_clin = profiles_with_clin[profiles_with_clin.n_bcells >= 1].dropna(subset=["clinical_progr"])
stat4, pval4 = mannwhitneyu(
    bcell_clin[bcell_clin.clinical_progr == 1]["tls_score"].values,
    bcell_clin[bcell_clin.clinical_progr == 0]["tls_score"].values,
    alternative="two-sided"
)
print(f"  TLS score vs clinical progression: p={pval4:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
for ax, (col, label, pval_val) in zip(axes, [
    ("psa_progr", "PSA recurrence", pval3),
    ("clinical_progr", "Clinical progression", pval4),
]):
    sub = profiles_with_clin[profiles_with_clin.n_bcells >= 1].dropna(subset=[col])
    sub[col] = sub[col].astype(str)
    sns.boxplot(data=sub, x=col, y="tls_score", ax=ax, palette="Set2")
    ax.set_title(f"TLS score vs {label}\np={pval_val:.3f}")
    ax.set_xlabel(label + " (0=no, 1=yes)")
    ax.set_ylabel("TLS score")

plt.tight_layout()
fig.savefig(OUT_FIG / "TLS_score_outcome.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 4: Treg:T ratio in B-cell neighborhoods ───────────────────────

# Compute Treg:T ratio per ROI
bcell_tc = profiles_with_clin[profiles_with_clin.n_bcells >= 1].copy()
# Avoid division by zero
tls_th = bcell_tc["tls_thelper_score"]
tls_tr = bcell_tc["tls_treg_score"]
tls_tc = bcell_tc["tls_score"]  # any T cell
bcell_tc["treg_to_t_ratio"] = tls_tr / (tls_tc + 1e-5)

# TLS maturity proxy: high B-T colocalization + low Treg ratio = potentially active TLS
bcell_tc["tls_activation_score"] = tls_th - tls_tr

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

ax = axes[0]
bc = bcell_tc.dropna(subset=["inflammation"])
sns.violinplot(data=bc, x="inflammation", y="treg_to_t_ratio", order=["no", "yes"],
               ax=ax, palette={"no": "#aec6cf", "yes": "#ff7f7f"}, inner="box")
ax.set_title("Treg:T ratio in B-cell neighborhoods\nby inflammation status")
ax.set_xlabel("Histological inflammation")
ax.set_ylabel("Treg/(any T cell) ratio")

ax = axes[1]
bc_gg = bcell_tc.dropna(subset=["gleason_grp"])
sns.boxplot(data=bc_gg, x="gleason_grp", y="tls_activation_score", ax=ax, palette="Blues")
ax.set_title("TLS activation score by Gleason grade group\n(helper T – Treg coloc with B cells)")
ax.set_xlabel("Gleason grade group")
ax.set_ylabel("T helper – Treg coloc score")

plt.tight_layout()
fig.savefig(OUT_FIG / "treg_ratio_and_tls_activation.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 5: B-cell density and TLS score scatter ────────────────────────

fig, ax = plt.subplots(figsize=(7, 5))
bc_all = profiles_with_clin.copy()
color_map = {"yes": "#ff7f7f", "no": "#aec6cf"}
for inf_val in ["yes", "no"]:
    sub = bc_all[bc_all.inflammation == inf_val].dropna(subset=["inflammation"])
    ax.scatter(sub["bcell_proportion"], sub["tls_score"], c=color_map.get(inf_val, "gray"),
               alpha=0.5, s=20, label=f"inflammation={inf_val}")
rho, rho_p = spearmanr(bc_all["bcell_proportion"].values, bc_all["tls_score"].values)
ax.set_xlabel("B-cell proportion")
ax.set_ylabel("TLS score")
ax.set_title(f"B-cell proportion vs TLS score\nSpearman rho={rho:.3f}, p={rho_p:.3e}")
ax.legend()
plt.tight_layout()
fig.savefig(OUT_FIG / "bcell_proportion_vs_tls_score.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Summary statistics ───────────────────────────────────────────────────

print("\n=== RQ3 Summary ===")
print(f"Total tumor ROIs analyzed: {len(profiles_with_clin)}")
print(f"ROIs with ≥1 B cell: {(profiles_with_clin.n_bcells >= 1).sum()}")
print(f"ROIs with ≥3 B cells: {(profiles_with_clin.n_bcells >= 3).sum()}")

print(f"\nTLS score (B cells with ≥1 T cell neighbor):")
print(f"  All ROIs with ≥1 B cell: mean={profiles_with_clin[profiles_with_clin.n_bcells>=1].tls_score.mean():.3f}")
infl_sub = profiles_with_clin[profiles_with_clin.n_bcells>=1].dropna(subset=["inflammation"])
print(f"  Inflamed: {infl_sub[infl_sub.inflammation=='yes'].tls_score.mean():.3f}")
print(f"  Non-inflamed: {infl_sub[infl_sub.inflammation=='no'].tls_score.mean():.3f}")
print(f"  MWU p-value (inflamed vs non-inflamed): {pval:.4f}")

print(f"\nPSA recurrence vs TLS score: p={pval3:.4f}")
print(f"Clinical progression vs TLS score: p={pval4:.4f}")

# Save summary
summary = profiles_with_clin[["n_bcells", "n_cells", "bcell_proportion",
                               "tls_score", "tls_thelper_score", "tls_treg_score",
                               "inflammation", "psa_progr", "clinical_progr",
                               "gleason_grp", "pat_id"]].copy()
summary.to_csv(OUT_TAB / "roi_tls_scores.csv")
print(f"\nOutputs written to {OUT_TAB}/ and {OUT_FIG}/")
