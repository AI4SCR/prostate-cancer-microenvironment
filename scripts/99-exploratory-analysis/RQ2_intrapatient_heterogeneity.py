"""
RQ2: Intra-patient Heterogeneity Between ROIs

For patients with multiple ROIs, quantifies:
1. Within-patient vs between-patient variance in cell composition (ICC)
2. Whether higher intra-patient heterogeneity tracks with clinical variables
3. Paired comparison of high-grade vs low-grade cores within the same patient
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu
from scipy.spatial.distance import jensenshannon
from statsmodels.stats.multitest import multipletests
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ2")
OUT_TAB = Path("output/tables/RQ2")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load data ─────────────────────────────────────────────────────────────

print("Loading data...")
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
clinical = pd.read_parquet(BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet")
tumor_ids = set(clinical.index[clinical.is_tumor == "yes"])

cell_meta = pd.concat(
    [pd.read_parquet(f, engine="fastparquet") for f in meta_files], ignore_index=False
)
cell_meta = cell_meta.loc[cell_meta.index.get_level_values("sample_id").isin(tumor_ids)]

# ROI-level meta_label proportions
meta_counts = cell_meta.groupby(["sample_id", "meta_label"]).size().unstack(fill_value=0)
meta_props = meta_counts.div(meta_counts.sum(axis=1), axis=0)

# Add clinical variables to ROI props
meta_props_clin = meta_props.join(
    clinical[["pat_id", "gleason_grp", "psa_progr", "clinical_progr", "psa_progr_time",
               "clinical_progr_time", "stromogenic_smc_loss_reactive_stroma_present"]]
)
print(f"  {len(meta_props_clin)} tumor ROIs")
print(f"  {meta_props_clin.pat_id.nunique()} patients")

# %% ─── Patient-level heterogeneity ──────────────────────────────────────────

print("Computing intra-patient heterogeneity...")

# CLR transformation
def clr(df: pd.DataFrame, pseudocount: float = 1e-5) -> pd.DataFrame:
    X = df.values + pseudocount
    log_X = np.log(X)
    log_gm = log_X.mean(axis=1, keepdims=True)
    return pd.DataFrame(log_X - log_gm, index=df.index, columns=df.columns)

ct_cols = [c for c in meta_props.columns]
props_clr = clr(meta_props)
props_clr_clin = props_clr.join(clinical[["pat_id", "gleason_grp", "psa_progr",
                                            "clinical_progr", "psa_progr_time",
                                            "clinical_progr_time"]])

# Within-patient SD per cell type
within_sd = props_clr_clin.groupby("pat_id")[ct_cols].std()
n_rois_per_pat = meta_props_clin.groupby("pat_id").size()
print(f"  Patients with ≥2 ROIs: {(n_rois_per_pat >= 2).sum()}")

# Aggregate: mean within-patient SD across all cell types = heterogeneity score
within_sd["mean_sd"] = within_sd.mean(axis=1)
within_sd = within_sd.join(n_rois_per_pat.rename("n_rois"))

# Join clinical
pat_clinical_sub = clinical.groupby("pat_id").first()[
    ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time", "gleason_grp"]
]
within_sd = within_sd.join(pat_clinical_sub)
within_sd.to_csv(OUT_TAB / "patient_intrapatient_heterogeneity.csv")

# %% ─── ICC per cell type ─────────────────────────────────────────────────────

# Intraclass correlation: ICC(1) = (MSB - MSW) / (MSB + (k-1)*MSW)
# where k = mean group size, MSB = mean square between patients, MSW = within

def compute_icc(props_clr_clin: pd.DataFrame, ct_col: str) -> float:
    """One-way ICC for a cell type across patients."""
    groups = props_clr_clin.dropna(subset=[ct_col, "pat_id"]).groupby("pat_id")[ct_col]
    group_data = [g.values for _, g in groups if len(g) >= 2]
    if len(group_data) < 5:
        return np.nan

    grand_mean = np.concatenate(group_data).mean()
    k_mean = np.mean([len(g) for g in group_data])
    n_groups = len(group_data)
    n_total = sum(len(g) for g in group_data)

    # SSB = sum over groups of n_j * (mean_j - grand_mean)^2
    ssb = sum(len(g) * (g.mean() - grand_mean)**2 for g in group_data)
    # SSW = sum over groups of sum (x_ij - mean_j)^2
    ssw = sum(((g - g.mean())**2).sum() for g in group_data)

    dfb = n_groups - 1
    dfw = n_total - n_groups
    if dfb <= 0 or dfw <= 0:
        return np.nan

    msb = ssb / dfb
    msw = ssw / dfw
    icc = (msb - msw) / (msb + (k_mean - 1) * msw) if (msb + (k_mean - 1) * msw) > 0 else np.nan
    return float(icc)

icc_scores = {}
for ct in ct_cols:
    icc_scores[ct] = compute_icc(props_clr_clin, ct)

icc_df = pd.Series(icc_scores, name="ICC").sort_values(ascending=False).to_frame()
icc_df.to_csv(OUT_TAB / "icc_per_cell_type.csv")

print("\nICC per cell type (sorted by ICC, descending = more stable within patients):")
print(icc_df.to_string())

# %% ─── JSD between ROIs within patients ─────────────────────────────────────

print("\nComputing pairwise JSD between ROIs within each patient...")

pat_jsd_records = []
for pat, group in meta_props_clin.groupby("pat_id"):
    if len(group) < 2:
        continue
    props_mat = group[ct_cols].values
    jsds = []
    for i in range(len(props_mat)):
        for j in range(i + 1, len(props_mat)):
            # JSD requires probability distributions
            p = props_mat[i] + 1e-10
            q = props_mat[j] + 1e-10
            p /= p.sum()
            q /= q.sum()
            jsds.append(jensenshannon(p, q))
    mean_jsd = np.mean(jsds) if jsds else np.nan
    pat_jsd_records.append({"pat_id": pat, "mean_jsd": mean_jsd, "n_rois": len(group)})

pat_jsd = pd.DataFrame(pat_jsd_records).set_index("pat_id")
pat_jsd = pat_jsd.join(pat_clinical_sub)
print(f"  {len(pat_jsd)} patients with ≥2 ROIs")

# %% ─── Figures ──────────────────────────────────────────────────────────────

print("Generating figures...")

# Figure 1: ICC per cell type (bar plot)
fig, ax = plt.subplots(figsize=(10, 6))
colors = ["steelblue" if v > 0.5 else ("orange" if v > 0.2 else "salmon") for v in icc_df.ICC]
ax.barh(icc_df.index, icc_df.ICC, color=colors)
ax.axvline(0, color="k", linewidth=0.5)
ax.axvline(0.5, color="steelblue", linewidth=0.8, linestyle="--", alpha=0.5, label="ICC=0.5")
ax.axvline(0.2, color="orange", linewidth=0.8, linestyle="--", alpha=0.5, label="ICC=0.2")
ax.set_xlabel("ICC (intraclass correlation coefficient)")
ax.set_title("Cell type stability within patients (ICC)\nHigh ICC = consistent across ROIs from same patient")
ax.legend(fontsize=9)
plt.tight_layout()
fig.savefig(OUT_FIG / "icc_per_cell_type.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# Figure 2: Mean within-patient SD by outcome
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
for ax, (col, label) in zip(axes, [
    ("psa_progr", "PSA recurrence"),
    ("clinical_progr", "Clinical progression"),
]):
    sub = within_sd[within_sd.n_rois >= 2].dropna(subset=["mean_sd", col])
    g0 = sub[sub[col] == 0]["mean_sd"].values
    g1 = sub[sub[col] == 1]["mean_sd"].values
    if len(g0) > 3 and len(g1) > 3:
        stat, pval = mannwhitneyu(g0, g1, alternative="two-sided")
    else:
        pval = np.nan
    sub_plot = sub.copy()
    sub_plot[col] = sub_plot[col].astype(str)
    sns.boxplot(data=sub_plot, x=col, y="mean_sd", ax=ax)
    ax.set_title(f"Intra-patient heterogeneity\nvs {label}\np={pval:.3f}" if not np.isnan(pval) else f"Intra-patient heterogeneity\nvs {label}")
    ax.set_xlabel(label + " (0=no, 1=yes)")
    ax.set_ylabel("Mean within-patient SD (CLR-transformed proportions)")

plt.tight_layout()
fig.savefig(OUT_FIG / "heterogeneity_vs_outcome.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# Figure 3: JSD by Gleason grade
fig, ax = plt.subplots(figsize=(8, 5))
gleason_map = {1.0: "GG1", 2.0: "GG2", 3.0: "GG3", 4.0: "GG4", 5.0: "GG5"}
sub = pat_jsd.dropna(subset=["mean_jsd", "gleason_grp"])
sub_plot = sub.copy()
sub_plot["gg"] = sub_plot["gleason_grp"].map(gleason_map)
rho, pval = spearmanr(sub["gleason_grp"].values, sub["mean_jsd"].values)
sns.boxplot(data=sub_plot, x="gg", y="mean_jsd", order=["GG1", "GG2", "GG3", "GG4", "GG5"], ax=ax)
ax.set_title(f"Intra-patient ROI divergence (Jensen-Shannon) by Gleason grade\nSpearman rho={rho:.3f}, p={pval:.3f}")
ax.set_xlabel("Gleason grade group (max per patient)")
ax.set_ylabel("Mean pairwise JSD between ROIs")
plt.tight_layout()
fig.savefig(OUT_FIG / "jsd_by_gleason.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# Figure 4: N ROIs per patient + ROI composition scatter
n_rois = n_rois_per_pat.value_counts().sort_index()
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
ax.bar(n_rois.index, n_rois.values)
ax.set_xlabel("Number of tumor ROIs per patient")
ax.set_ylabel("Number of patients")
ax.set_title("ROI sampling per patient")

# PCA of composition space
from sklearn.decomposition import PCA
ax = axes[1]
pca = PCA(n_components=2)
props_clr_arr = props_clr.values
pca_coords = pca.fit_transform(props_clr_arr)
pca_df = pd.DataFrame(pca_coords, index=props_clr.index, columns=["PC1", "PC2"])
pca_df = pca_df.join(clinical[["pat_id", "gleason_grp"]])
# Color by Gleason
scatter = ax.scatter(pca_df.PC1, pca_df.PC2,
                      c=pca_df.gleason_grp, cmap="RdYlBu_r",
                      alpha=0.5, s=15)
plt.colorbar(scatter, ax=ax, label="Gleason grade group")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA of CLR-transformed cell composition\ncolored by Gleason grade group")

plt.tight_layout()
fig.savefig(OUT_FIG / "roi_pca_composition.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Summary ───────────────────────────────────────────────────────────────

print("\n=== RQ2 Summary ===")
print(f"Patients with ≥2 tumor ROIs: {(n_rois_per_pat >= 2).sum()}")
print(f"\nTop ICC cell types (most stable within patients):")
print(icc_df.head(5).to_string())
print(f"\nBottom ICC cell types (most variable within patients):")
print(icc_df.tail(5).to_string())

# JSD vs Gleason
sub = pat_jsd.dropna(subset=["mean_jsd", "gleason_grp"])
rho, pval = spearmanr(sub["gleason_grp"].values, sub["mean_jsd"].values)
print(f"\nIntra-patient JSD vs Gleason: rho={rho:.3f}, p={pval:.4f}")

print(f"\nOutputs saved to {OUT_TAB}/ and {OUT_FIG}/")
