"""
RQ8: Paired High-Grade vs Low-Grade TMA Core Analysis Within Patients

Uses the gleason_pattern_tma_core annotation ("high"/"low") to compare TME composition
and spatial diversity between matched high- and low-grade cores from the same patient.
Wilcoxon signed-rank paired tests control for patient-level confounders.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ8")
OUT_TAB = Path("output/tables/RQ8")
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

# CLR transformation
def clr(df: pd.DataFrame, pseudocount: float = 1e-5) -> pd.DataFrame:
    X = df.values + pseudocount
    log_X = np.log(X)
    return pd.DataFrame(log_X - log_X.mean(axis=1, keepdims=True), index=df.index, columns=df.columns)

props_clr = clr(meta_props)

# Attach clinical variables including gleason_pattern_tma_core
roi_info = clinical[["pat_id", "gleason_grp", "gleason_pattern_tma_core",
                       "psa_progr", "clinical_progr"]]
props_clin = props_clr.join(roi_info)

print(f"  {len(props_clin)} tumor ROIs")
print(f"  gleason_pattern_tma_core distribution:")
print(props_clin.gleason_pattern_tma_core.value_counts().to_dict())

# %% ─── Identify paired patients ──────────────────────────────────────────────

# Keep only "high" and "low" labeled ROIs
props_hl = props_clin[props_clin.gleason_pattern_tma_core.isin(["high", "low"])].copy()
ct_cols = [c for c in meta_props.columns]

# Patients with at least 1 high and 1 low ROI
pat_counts = props_hl.groupby(["pat_id", "gleason_pattern_tma_core"]).size().unstack(fill_value=0)
paired_pats = pat_counts[(pat_counts.get("high", 0) >= 1) & (pat_counts.get("low", 0) >= 1)].index
print(f"\nPatients with ≥1 high AND ≥1 low grade core: {len(paired_pats)}")

# Average within each patient×grade category
paired_props = (
    props_hl[props_hl.pat_id.isin(paired_pats)]
    .groupby(["pat_id", "gleason_pattern_tma_core"])[ct_cols]
    .mean()
)

# Reshape to (pat_id × [high_ctX, low_ctX])
high_props = paired_props.xs("high", level=1).add_prefix("high_")
low_props  = paired_props.xs("low",  level=1).add_prefix("low_")
paired_df = high_props.join(low_props)

# Also load Shannon diversity scores (from RQ5)
shannon_file = Path("output/tables/RQ5/roi_shannon_diversity.csv")
if shannon_file.exists():
    shannon_df = pd.read_csv(shannon_file, index_col=0)
    shannon_df["pat_id"] = clinical.loc[shannon_df.index.intersection(clinical.index), "pat_id"]
    shannon_df["grade"] = clinical.loc[shannon_df.index.intersection(clinical.index), "gleason_pattern_tma_core"]
    shannon_hl = shannon_df[shannon_df.grade.isin(["high", "low"])].copy()
    shannon_hl = shannon_hl[shannon_hl.pat_id.isin(paired_pats)]
    shannon_paired = shannon_hl.groupby(["pat_id", "grade"])[["h_global", "h_local_mean", "h_local_sd"]].mean()
    sha_high = shannon_paired.xs("high", level=1).add_prefix("sha_high_")
    sha_low  = shannon_paired.xs("low",  level=1).add_prefix("sha_low_")
    paired_df = paired_df.join(sha_high).join(sha_low)
    HAS_SHANNON = True
else:
    HAS_SHANNON = False
    print("  Shannon file not found; skipping spatial diversity comparison")

print(f"Paired DataFrame: {len(paired_df)} patients × {paired_df.shape[1]} features")
paired_df.to_csv(OUT_TAB / "paired_high_low_props.csv")

# %% ─── Paired Wilcoxon signed-rank tests ────────────────────────────────────

print("\nPaired Wilcoxon tests (high vs low grade cores, within-patient)...")

wilcoxon_results = []
for ct in ct_cols:
    high_col = f"high_{ct}"
    low_col  = f"low_{ct}"
    if high_col not in paired_df.columns or low_col not in paired_df.columns:
        continue
    diff = paired_df[high_col].values - paired_df[low_col].values
    # Drop NaN
    diff = diff[~np.isnan(diff)]
    if len(diff) < 5 or (diff == 0).all():
        continue
    try:
        stat, pval = wilcoxon(diff, alternative="two-sided")
    except ValueError:
        continue
    mean_diff = diff.mean()
    # Rank-biserial correlation as effect size for Wilcoxon
    n = len(diff)
    rbc = 1 - 2 * stat / (n * (n + 1) / 2) if n > 0 else np.nan
    wilcoxon_results.append({
        "cell_type": ct, "mean_diff_high_minus_low": mean_diff,
        "wilcoxon_stat": stat, "pval": pval, "rbc": rbc, "n_pairs": n,
    })

wx_df = pd.DataFrame(wilcoxon_results)
wx_df["fdr"] = multipletests(wx_df["pval"].fillna(1), method="fdr_bh")[1]
wx_df = wx_df.sort_values("fdr")
wx_df.to_csv(OUT_TAB / "paired_wilcoxon_results.csv", index=False)

print(f"Significant (FDR < 0.1): {(wx_df.fdr < 0.1).sum()}")
print(wx_df[wx_df.fdr < 0.1][["cell_type","mean_diff_high_minus_low","rbc","fdr"]].to_string(index=False))

# %% ─── Figure 1: Lollipop plot — high minus low ──────────────────────────────

fig, ax = plt.subplots(figsize=(10, 7))
sig_mask = wx_df["fdr"] < 0.1
colors = ["red" if s else "gray" for s in sig_mask]

ax.hlines(wx_df["cell_type"], 0, wx_df["mean_diff_high_minus_low"],
          color=colors, linewidth=1.5, alpha=0.8)
ax.scatter(wx_df["mean_diff_high_minus_low"], wx_df["cell_type"],
           color=colors, s=50, zorder=5)
ax.axvline(0, color="black", linewidth=0.8, linestyle="--")

# Annotate significant ones
for _, row in wx_df[sig_mask].iterrows():
    ax.annotate(f"FDR={row['fdr']:.3f}",
                (row["mean_diff_high_minus_low"], row["cell_type"]),
                fontsize=7, xytext=(5, 0), textcoords="offset points", va="center")

ax.set_xlabel("Mean CLR proportion (high grade – low grade)")
ax.set_title(f"Paired comparison: high vs low Gleason TMA cores within patient\n"
             f"n={len(paired_df)} patients; red = FDR < 0.1")
ax.set_ylabel("Cell type (meta-label)")
plt.tight_layout()
fig.savefig(OUT_FIG / "paired_lollipop.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 2: Slope plot for top significant cell types ──────────────────

top_ct = wx_df[wx_df.fdr < 0.1]["cell_type"].tolist()

if top_ct:
    n_ct = len(top_ct)
    ncols = min(3, n_ct)
    nrows = (n_ct + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes_flat = np.array(axes).ravel() if n_ct > 1 else [axes]

    for ax, ct in zip(axes_flat, top_ct):
        high_col = f"high_{ct}"
        low_col  = f"low_{ct}"
        sub = paired_df[[high_col, low_col]].dropna()
        fdr_val = wx_df[wx_df.cell_type == ct]["fdr"].values[0]
        mean_diff = wx_df[wx_df.cell_type == ct]["mean_diff_high_minus_low"].values[0]

        for _, row in sub.iterrows():
            color = "red" if row[high_col] > row[low_col] else "steelblue"
            ax.plot([0, 1], [row[low_col], row[high_col]], color=color, alpha=0.3, linewidth=0.8)

        # Mean line
        ax.plot([0, 1], [sub[low_col].mean(), sub[high_col].mean()],
                color="black", linewidth=2.5, linestyle="--", label="Mean")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Low grade", "High grade"])
        ax.set_ylabel("CLR proportion")
        ax.set_title(f"{ct}\nFDR={fdr_val:.3f}, Δ={mean_diff:.3f}", fontsize=8)

    for ax in axes_flat[n_ct:]:
        ax.set_visible(False)

    plt.suptitle("Paired high vs low grade core: CLR cell type proportions", fontsize=11)
    plt.tight_layout()
    fig.savefig(OUT_FIG / "paired_slope_plot.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

# %% ─── Figure 3: Spatial diversity in high vs low grade cores ───────────────

if HAS_SHANNON:
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    for ax, (sha_metric, label) in zip(axes, [
        ("h_global", "Global Shannon"),
        ("h_local_mean", "Mean local Shannon"),
        ("h_local_sd", "SD local Shannon"),
    ]):
        hc = f"sha_high_{sha_metric}"
        lc = f"sha_low_{sha_metric}"
        if hc not in paired_df.columns or lc not in paired_df.columns:
            continue
        sub = paired_df[[hc, lc]].dropna()
        diff = sub[hc].values - sub[lc].values
        try:
            stat, pval = wilcoxon(diff, alternative="two-sided")
        except ValueError:
            pval = np.nan

        melt = pd.DataFrame({
            "grade": ["Low grade"] * len(sub) + ["High grade"] * len(sub),
            label: list(sub[lc]) + list(sub[hc]),
        })
        sns.boxplot(data=melt, x="grade", y=label, order=["Low grade", "High grade"], ax=ax,
                    palette={"Low grade": "#aec6cf", "High grade": "#ff7f7f"})
        ax.set_title(f"{label}\npaired Wilcoxon p={pval:.4f}" if not np.isnan(pval) else label)
        ax.set_xlabel("")

    plt.suptitle("Spatial diversity: high vs low grade cores within patients", fontsize=11)
    plt.tight_layout()
    fig.savefig(OUT_FIG / "paired_shannon_diversity.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

# %% ─── Figure 4: Volcano plot ────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(8, 6))
x = wx_df["rbc"].values
y = -np.log10(wx_df["fdr"].values + 1e-10)
colors = ["red" if f < 0.1 else "gray" for f in wx_df["fdr"]]
ax.scatter(x, y, c=colors, s=40, alpha=0.8)
ax.axhline(-np.log10(0.1), color="red", linestyle="--", linewidth=0.8, alpha=0.6)
ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)

for _, row in wx_df[wx_df.fdr < 0.2].iterrows():
    ax.annotate(row["cell_type"], (row["rbc"], -np.log10(row["fdr"] + 1e-10)),
                fontsize=7, xytext=(5, 2), textcoords="offset points")

ax.set_xlabel("Rank-biserial correlation (high > low if positive)")
ax.set_ylabel("-log10(FDR)")
ax.set_title(f"Paired high vs low grade volcano\n(n={len(paired_df)} patients)")
plt.tight_layout()
fig.savefig(OUT_FIG / "paired_volcano.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Summary ───────────────────────────────────────────────────────────────

print(f"\n=== RQ8 Summary ===")
print(f"Patients with paired high/low cores: {len(paired_pats)}")
print(f"Patients in analysis after averaging: {len(paired_df)}")
print(f"\nTop paired differences (high minus low grade):")
print(wx_df[["cell_type","mean_diff_high_minus_low","rbc","fdr"]].head(10).to_string(index=False))

if HAS_SHANNON:
    for sha_metric in ["h_global", "h_local_mean", "h_local_sd"]:
        hc = f"sha_high_{sha_metric}"
        lc = f"sha_low_{sha_metric}"
        if hc in paired_df.columns and lc in paired_df.columns:
            sub = paired_df[[hc, lc]].dropna()
            diff = sub[hc].values - sub[lc].values
            try:
                stat, pval = wilcoxon(diff)
                print(f"\n{sha_metric}: mean diff={diff.mean():.4f}, Wilcoxon p={pval:.4f}")
            except ValueError:
                pass

print(f"\nOutputs → {OUT_TAB}/ and {OUT_FIG}/")
