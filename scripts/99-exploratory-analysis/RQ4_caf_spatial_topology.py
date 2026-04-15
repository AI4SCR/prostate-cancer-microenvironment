"""
RQ4: CAF Subtype Spatial Topology Relative to Tumor Glands

Computes proximity of CAF1 (CD105+) and CAF2 (AR+) subtypes to luminal epithelium,
tests permutation-based interaction enrichment, and associates CAF–epithelial proximity
scores with clinical outcome.
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
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ4")
OUT_TAB = Path("output/tables/RQ4")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load data ─────────────────────────────────────────────────────────────

print("Loading data...")
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
spatial_files = sorted((BASE_DIR / "02_processed/features/spatial/filtered-annotated").glob("*.parquet"))

meta_by_id = {f.stem: f for f in meta_files}
spatial_by_id = {f.stem: f for f in spatial_files}

clinical = pd.read_parquet(BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet")
tumor_ids = set(clinical.index[clinical.is_tumor == "yes"])
sample_ids = sorted(set(meta_by_id) & set(spatial_by_id) & tumor_ids)
print(f"  {len(sample_ids)} tumor ROIs")

# %% ─── Cell type definitions ─────────────────────────────────────────────────

CAF1_PLUS_LABEL = "stromal-CAF1(CD105+)"
CAF1_MINUS_LABEL = "stromal-CAF1(CD105-)"
CAF1_EGR1_LABEL = "stromal-CAF1(CD105-EGR1+)"
CAF2_AR_PLUS_LABEL = "stromal-CAF2(AR+)"
CAF2_AR_CES1_LABEL = "stromal-CAF2(AR+CES1+)"
CAF2_AR_EGR1_LABEL = "stromal-CAF2(AR+EGR1+)"
CAF2_AR_MINUS_LABEL = "stromal-CAF2(AR-)"

CAF1_LABELS = {CAF1_PLUS_LABEL, CAF1_MINUS_LABEL, CAF1_EGR1_LABEL}
CAF2_LABELS = {CAF2_AR_PLUS_LABEL, CAF2_AR_CES1_LABEL, CAF2_AR_EGR1_LABEL, CAF2_AR_MINUS_LABEL}
CAF2_AR_POS_LABELS = {CAF2_AR_PLUS_LABEL, CAF2_AR_CES1_LABEL, CAF2_AR_EGR1_LABEL}

LUMINAL_LABELS = {
    "epithelial-luminal", "epithelial-luminal(ERG+)", "epithelial-luminal(p53+)",
    "epithelial-luminal(Ki67+)", "epithelial-(ERG+CD44+)"
}
TUMOR_LABELS = LUMINAL_LABELS | {"epithelial-transient"}

# %% ─── Per-ROI proximity analysis ───────────────────────────────────────────

print("Computing CAF proximity to luminal epithelium per ROI...")


def compute_caf_proximity(spatial_df: pd.DataFrame, meta_df: pd.DataFrame) -> dict:
    """
    Compute mean nearest-neighbor distances from CAF subtypes to luminal epithelium.
    Also compute interaction enrichment within radius-32 neighborhood.
    """
    common_obj = spatial_df.index.intersection(meta_df.index)
    if len(common_obj) < 5:
        return {}

    spatial_df = spatial_df.loc[common_obj]
    meta_df = meta_df.loc[common_obj]
    coords = spatial_df[["x", "y"]].values
    labels = meta_df["label"].values

    luminal_mask = np.array([l in LUMINAL_LABELS for l in labels])
    caf1_plus_mask = labels == CAF1_PLUS_LABEL
    caf2_ar_mask = np.array([l in CAF2_AR_POS_LABELS for l in labels])
    caf2_ar_minus_mask = labels == CAF2_AR_MINUS_LABEL
    all_caf1_mask = np.array([l in CAF1_LABELS for l in labels])
    all_caf2_mask = np.array([l in CAF2_LABELS for l in labels])

    result = {
        "n_cells": len(labels),
        "n_luminal": luminal_mask.sum(),
        "n_caf1_plus": caf1_plus_mask.sum(),
        "n_caf2_ar": caf2_ar_mask.sum(),
        "n_caf2_ar_minus": caf2_ar_minus_mask.sum(),
        "n_caf1_all": all_caf1_mask.sum(),
        "n_caf2_all": all_caf2_mask.sum(),
    }

    # Proportions
    n = len(labels)
    result["prop_caf1_plus"] = caf1_plus_mask.sum() / n
    result["prop_caf2_ar"] = caf2_ar_mask.sum() / n
    result["prop_caf2_ar_minus"] = caf2_ar_minus_mask.sum() / n
    result["prop_luminal"] = luminal_mask.sum() / n

    if luminal_mask.sum() == 0:
        return result

    luminal_coords = coords[luminal_mask]
    luminal_tree = cKDTree(luminal_coords)

    RADIUS = 32  # pixels

    for caf_name, caf_mask in [
        ("caf1_plus", caf1_plus_mask),
        ("caf2_ar", caf2_ar_mask),
        ("caf2_ar_minus", caf2_ar_minus_mask),
        ("caf1_all", all_caf1_mask),
        ("caf2_all", all_caf2_mask),
    ]:
        if caf_mask.sum() == 0:
            continue
        caf_coords = coords[caf_mask]

        # Mean nearest-neighbor distance to luminal
        nn_dists, _ = luminal_tree.query(caf_coords, k=1)
        result[f"{caf_name}_mean_nn_dist_to_luminal"] = nn_dists.mean()
        result[f"{caf_name}_median_nn_dist_to_luminal"] = np.median(nn_dists)

        # Fraction of CAF cells within radius of any luminal cell (periglandular index)
        counts_in_radius = luminal_tree.query_ball_point(caf_coords, r=RADIUS)
        n_adjacent = sum(len(c) > 0 for c in counts_in_radius)
        result[f"{caf_name}_periglandular_fraction"] = n_adjacent / caf_mask.sum()

    # Interaction enrichment: observed vs expected (shuffled)
    # Count CAF1+ - luminal contacts in radius-32 neighborhood
    all_tree = cKDTree(coords)

    def count_cross_contacts(mask_a, mask_b, radius=RADIUS):
        if mask_a.sum() == 0 or mask_b.sum() == 0:
            return 0, 0
        coords_a = coords[mask_a]
        tree_b = cKDTree(coords[mask_b])
        counts = tree_b.query_ball_point(coords_a, r=radius)
        return sum(len(c) for c in counts), mask_a.sum()

    obs_caf1plus_lum, n_caf1plus = count_cross_contacts(caf1_plus_mask, luminal_mask)
    obs_caf2ar_lum, n_caf2ar = count_cross_contacts(caf2_ar_mask, luminal_mask)

    result["obs_caf1plus_luminal_contacts"] = obs_caf1plus_lum
    result["obs_caf2ar_luminal_contacts"] = obs_caf2ar_lum

    # Expected: shuffle labels (5 permutations for speed)
    n_perm = 20
    exp_caf1plus_lum = []
    exp_caf2ar_lum = []
    perm_labels = labels.copy()
    for _ in range(n_perm):
        np.random.shuffle(perm_labels)
        perm_caf1p = perm_labels == CAF1_PLUS_LABEL
        perm_caf2ar = np.array([l in CAF2_AR_POS_LABELS for l in perm_labels])
        perm_lum = np.array([l in LUMINAL_LABELS for l in perm_labels])
        c1, _ = count_cross_contacts(perm_caf1p, perm_lum)
        c2, _ = count_cross_contacts(perm_caf2ar, perm_lum)
        exp_caf1plus_lum.append(c1)
        exp_caf2ar_lum.append(c2)

    exp_c1 = np.mean(exp_caf1plus_lum) if exp_caf1plus_lum else 1
    exp_c2 = np.mean(exp_caf2ar_lum) if exp_caf2ar_lum else 1
    result["enrichment_caf1plus_luminal"] = obs_caf1plus_lum / max(exp_c1, 1e-5)
    result["enrichment_caf2ar_luminal"] = obs_caf2ar_lum / max(exp_c2, 1e-5)

    return result


np.random.seed(42)
roi_caf_profiles = []
for i, sid in enumerate(sample_ids):
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(sample_ids)}")
    spatial_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    spatial_df = spatial_df.loc[sid] if sid in spatial_df.index.get_level_values("sample_id") else spatial_df
    meta_df = meta_df.loc[sid] if sid in meta_df.index.get_level_values("sample_id") else meta_df
    profile = compute_caf_proximity(spatial_df, meta_df)
    profile["sample_id"] = sid
    roi_caf_profiles.append(profile)

caf_df = pd.DataFrame(roi_caf_profiles).set_index("sample_id")
print(f"  {len(caf_df)} ROIs processed")

caf_df.to_csv(OUT_TAB / "roi_caf_proximity_profiles.csv")

# Merge with clinical
caf_clin = caf_df.join(clinical[["gleason_grp", "psa_progr", "psa_progr_time",
                                    "clinical_progr", "clinical_progr_time",
                                    "stromogenic_smc_loss_reactive_stroma_present",
                                    "cribriform", "inflammation", "pat_id"]], how="left")

# %% ─── Figure 1: CAF1+ vs CAF2-AR+ proximity to luminal epithelium ──────────

print("Generating proximity figures...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 1a. Mean NN distance: CAF1+ vs CAF2-AR+
both = caf_clin.dropna(subset=["caf1_plus_mean_nn_dist_to_luminal", "caf2_ar_mean_nn_dist_to_luminal"])
both_plot = both[both.n_caf1_plus >= 3][both.n_caf2_ar >= 3]

ax = axes[0]
ax.scatter(both_plot["caf1_plus_mean_nn_dist_to_luminal"],
           both_plot["caf2_ar_mean_nn_dist_to_luminal"],
           alpha=0.4, s=20)
ax.plot([0, 200], [0, 200], "k--", linewidth=0.8)
ax.set_xlabel("CAF1 CD105+ mean NN dist to luminal epithelium (px)")
ax.set_ylabel("CAF2 AR+ mean NN dist to luminal epithelium (px)")
ax.set_title("CAF proximity to luminal epithelium\n(one point per ROI)")

# 1b. Periglandular fraction comparison
melt_cols = {
    "caf1_plus_periglandular_fraction": "CAF1-CD105+",
    "caf2_ar_periglandular_fraction": "CAF2-AR+",
    "caf2_ar_minus_periglandular_fraction": "CAF2-AR-",
}
plot_data = []
for col, label in melt_cols.items():
    vals = caf_clin.dropna(subset=[col])
    for v in vals[col].values:
        plot_data.append({"CAF type": label, "Periglandular fraction": v})

pf_df = pd.DataFrame(plot_data)
ax = axes[1]
sns.boxplot(data=pf_df, x="CAF type", y="Periglandular fraction", ax=ax)
ax.set_title("Fraction of CAF cells within 32px\nof luminal epithelium (periglandular index)")
ax.set_ylabel("Periglandular fraction")
ax.tick_params(axis="x", rotation=15)

# 1c. Enrichment scores
enrich_plot = caf_clin.dropna(subset=["enrichment_caf1plus_luminal", "enrichment_caf2ar_luminal"])
ax = axes[2]
ax.hist(enrich_plot["enrichment_caf1plus_luminal"].clip(0, 5), bins=30, alpha=0.6, label="CAF1-CD105+", color="steelblue")
ax.hist(enrich_plot["enrichment_caf2ar_luminal"].clip(0, 5), bins=30, alpha=0.6, label="CAF2-AR+", color="salmon")
ax.axvline(1.0, color="k", linestyle="--", linewidth=0.8)
ax.set_xlabel("CAF–Luminal contact enrichment\n(observed/expected)")
ax.set_ylabel("N ROIs")
ax.set_title("Permutation-based contact enrichment\nCAF1+ and CAF2-AR+ vs luminal epithelium")
ax.legend()

plt.tight_layout()
fig.savefig(OUT_FIG / "caf_proximity_overview.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 2: CAF proximity by Gleason grade ─────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
gleason_map = {1.0: "GG1", 2.0: "GG2", 3.0: "GG3", 4.0: "GG4", 5.0: "GG5"}

for ax, (col, label) in zip(axes, [
    ("caf2_ar_mean_nn_dist_to_luminal", "CAF2-AR+ mean NN dist to luminal (px)"),
    ("caf2_ar_periglandular_fraction", "CAF2-AR+ periglandular fraction"),
]):
    sub = caf_clin.dropna(subset=[col, "gleason_grp"])
    sub = sub.copy()
    sub["gg_label"] = sub["gleason_grp"].map(gleason_map)
    sns.boxplot(data=sub, x="gg_label", y=col, ax=ax, palette="Blues",
                order=["GG1", "GG2", "GG3", "GG4", "GG5"])
    rho, pval = spearmanr(sub["gleason_grp"].values, sub[col].values)
    ax.set_title(f"{label}\nSpearman rho={rho:.3f}, p={pval:.3e}")
    ax.set_xlabel("Gleason grade group")
    ax.set_ylabel(label)

plt.tight_layout()
fig.savefig(OUT_FIG / "caf2_proximity_by_gleason.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 3: CAF proximity by stromogenic and cribriform ─────────────────

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for row_idx, (col, col_label) in enumerate([
    ("caf2_ar_mean_nn_dist_to_luminal", "CAF2-AR+ mean NN dist to luminal (px)"),
    ("caf1_plus_periglandular_fraction", "CAF1-CD105+ periglandular fraction"),
]):
    for col_idx, (clin_col, clin_label) in enumerate([
        ("stromogenic_smc_loss_reactive_stroma_present", "Stromogenic"),
        ("cribriform", "Cribriform"),
    ]):
        ax = axes[row_idx, col_idx]
        sub = caf_clin.dropna(subset=[col, clin_col])
        g0 = sub[sub[clin_col] == "no"][col].values
        g1 = sub[sub[clin_col] == "yes"][col].values
        if len(g0) > 3 and len(g1) > 3:
            stat, pval = mannwhitneyu(g0, g1, alternative="two-sided")
            n0, n1 = len(g0), len(g1)
            rbc = 1 - (2 * stat) / (n0 * n1)
        else:
            pval, rbc = np.nan, np.nan
        sns.boxplot(data=sub, x=clin_col, y=col, ax=ax, order=["no", "yes"],
                    palette={"no": "#aec6cf", "yes": "#ff7f7f"})
        ax.set_title(f"{col_label} by {clin_label}\np={pval:.3f}, rbc={rbc:.3f}" if not np.isnan(pval) else f"{col_label} by {clin_label}")
        ax.set_xlabel(clin_label)
        ax.set_ylabel(col_label)

plt.tight_layout()
fig.savefig(OUT_FIG / "caf_proximity_stromogenic_cribriform.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 4: CAF1+ periglandular index vs PSA recurrence ─────────────────

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

for ax, (col, clin_col, clin_label) in zip(axes, [
    ("caf1_plus_periglandular_fraction", "psa_progr", "PSA recurrence"),
    ("enrichment_caf1plus_luminal", "psa_progr", "PSA recurrence"),
]):
    sub = caf_clin.dropna(subset=[col, clin_col])
    g0 = sub[sub[clin_col] == 0][col].values
    g1 = sub[sub[clin_col] == 1][col].values
    if len(g0) > 3 and len(g1) > 3:
        stat, pval = mannwhitneyu(g0, g1, alternative="two-sided")
        n0, n1 = len(g0), len(g1)
        rbc = 1 - (2 * stat) / (n0 * n1)
    else:
        pval, rbc = np.nan, np.nan
    sub = sub.copy()
    sub[clin_col] = sub[clin_col].astype(str)
    sns.boxplot(data=sub, x=clin_col, y=col, ax=ax, palette="Set2")
    ax.set_title(f"{col} vs {clin_label}\np={pval:.3f}, rbc={rbc:.3f}" if not np.isnan(pval) else f"{col} vs {clin_label}")
    ax.set_xlabel(clin_label + " (0=no, 1=yes)")
    ax.set_ylabel(col)

plt.tight_layout()
fig.savefig(OUT_FIG / "caf1_periglandular_vs_outcome.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Statistical summary ───────────────────────────────────────────────────

print("\n=== RQ4 Statistical Summary ===")

# Compare CAF1+ vs CAF2-AR+ proximity
paired = caf_clin.dropna(subset=["caf1_plus_mean_nn_dist_to_luminal", "caf2_ar_mean_nn_dist_to_luminal"])
paired = paired[["caf1_plus_mean_nn_dist_to_luminal", "caf2_ar_mean_nn_dist_to_luminal"]].dropna()
print(f"\nMean NN distance to luminal epithelium (ROIs with both CAF types):")
print(f"  CAF1-CD105+: {paired.caf1_plus_mean_nn_dist_to_luminal.mean():.1f} px (median {paired.caf1_plus_mean_nn_dist_to_luminal.median():.1f})")
print(f"  CAF2-AR+:    {paired.caf2_ar_mean_nn_dist_to_luminal.mean():.1f} px (median {paired.caf2_ar_mean_nn_dist_to_luminal.median():.1f})")
if len(paired) > 5:
    stat, pval = wilcoxon(paired.caf1_plus_mean_nn_dist_to_luminal, paired.caf2_ar_mean_nn_dist_to_luminal)
    print(f"  Wilcoxon signed-rank p={pval:.4e} (CAF1+ vs CAF2-AR+ NN distance, n={len(paired)} ROIs)")

print("\nPeriglandular fraction:")
print(f"  CAF1-CD105+: {caf_clin.caf1_plus_periglandular_fraction.mean():.3f}")
print(f"  CAF2-AR+:    {caf_clin.caf2_ar_periglandular_fraction.mean():.3f}")

# Gleason correlation for CAF2-AR+ proximity
sub = caf_clin.dropna(subset=["caf2_ar_mean_nn_dist_to_luminal", "gleason_grp"])
rho, pval = spearmanr(sub.gleason_grp.values, sub.caf2_ar_mean_nn_dist_to_luminal.values)
print(f"\nCAF2-AR+ NN dist vs Gleason: Spearman rho={rho:.3f}, p={pval:.3e}")

sub2 = caf_clin.dropna(subset=["caf2_ar_periglandular_fraction", "gleason_grp"])
rho2, pval2 = spearmanr(sub2.gleason_grp.values, sub2.caf2_ar_periglandular_fraction.values)
print(f"CAF2-AR+ periglandular fraction vs Gleason: Spearman rho={rho2:.3f}, p={pval2:.3e}")

# Enrichment comparison
enr = caf_clin.dropna(subset=["enrichment_caf1plus_luminal", "enrichment_caf2ar_luminal"])
print(f"\nContact enrichment (observed/expected):")
print(f"  CAF1-CD105+: mean={enr.enrichment_caf1plus_luminal.mean():.2f}, median={enr.enrichment_caf1plus_luminal.median():.2f}")
print(f"  CAF2-AR+:    mean={enr.enrichment_caf2ar_luminal.mean():.2f}, median={enr.enrichment_caf2ar_luminal.median():.2f}")
stat, pval = wilcoxon(enr.enrichment_caf1plus_luminal.clip(0, 20), enr.enrichment_caf2ar_luminal.clip(0, 20))
print(f"  Wilcoxon signed-rank p={pval:.4e} (enrichment difference, n={len(enr)})")

print(f"\nOutputs saved to {OUT_TAB}/ and {OUT_FIG}/")
