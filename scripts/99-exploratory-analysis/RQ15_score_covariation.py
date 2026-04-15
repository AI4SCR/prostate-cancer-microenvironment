"""
RQ15: Co-Variation Between Patient-Level Spatial Scores

Builds a correlation matrix of all patient-level spatial TME scores and identifies which
properties co-vary. Visualizes the patient TME landscape and tests whether co-varying
score combinations better predict survival than individual scores.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from statsmodels.stats.multitest import multipletests
from lifelines import CoxPHFitter
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ15")
OUT_TAB = Path("output/tables/RQ15")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Collect all patient-level scores ─────────────────────────────────────

print("Collecting patient-level scores from prior analyses...")
clinical = pd.read_parquet(
    BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet"
)

def pat_id_map(index: pd.Index) -> pd.Series:
    return clinical.loc[index.intersection(clinical.index), "pat_id"]


scores: dict[str, pd.Series] = {}

# ── RQ5: Shannon diversity ──
rq5 = pd.read_csv("output/tables/RQ5/roi_shannon_diversity.csv", index_col=0)
rq5.index = rq5.index.astype(str)
rq5["pat_id"] = pat_id_map(rq5.index)
scores["h_global"]     = rq5.groupby("pat_id")["h_global"].mean()
scores["h_local_mean"] = rq5.groupby("pat_id")["h_local_mean"].mean()

# ── RQ4: CAF proximity ──
rq4 = pd.read_csv("output/tables/RQ4/roi_caf_proximity_profiles.csv", index_col=0)
rq4.index = rq4.index.astype(str)
rq4["pat_id"] = pat_id_map(rq4.index)
scores["caf1_periglandular"] = rq4.groupby("pat_id")["caf1_plus_periglandular_fraction"].mean()
scores["caf2_ar_periglandular"] = rq4.groupby("pat_id")["caf2_ar_periglandular_fraction"].mean()
scores["caf1_nn_dist"]  = rq4.groupby("pat_id")["caf1_plus_mean_nn_dist_to_luminal"].mean()

# ── RQ10: B-cell clusters ──
rq10 = pd.read_csv("output/tables/RQ10/roi_bcell_cluster_metrics.csv", index_col=0)
rq10.index = rq10.index.astype(str)
rq10["pat_id"] = pat_id_map(rq10.index)
scores["bcell_max_cluster"]      = rq10.groupby("pat_id")["max_cluster_size"].max()
scores["bcell_bb_density"]       = rq10.groupby("pat_id")["bb_contact_density"].mean()
scores["bcell_frac_large"]       = rq10.groupby("pat_id")["frac_in_large_cluster"].mean()

# ── RQ9: Niche usage (top 6 most variable niches) ──
rq9 = pd.read_csv("output/tables/RQ9/roi_niche_usage.csv", index_col=0)
rq9.index = rq9.index.astype(str)
rq9["pat_id"] = pat_id_map(rq9.index)
niche_cols = [c for c in rq9.columns if c.startswith("niche_")]
pat_niche = rq9.groupby("pat_id")[niche_cols].mean()
# Keep niches with highest between-patient variance
top_niche_cols = pat_niche.var().nlargest(8).index.tolist()
for c in top_niche_cols:
    scores[c] = pat_niche[c]

# ── RQ12: PC scores ──
pc_path = Path("output/tables/RQ12/patient_pc_scores.csv")
if pc_path.exists():
    pc_df = pd.read_csv(pc_path, index_col=0)
    for pc in ["PC1", "PC2", "PC3"]:
        if pc in pc_df.columns:
            scores[pc] = pc_df[pc]
    print("  RQ12 PC scores loaded")
else:
    print("  WARNING: RQ12 not yet run — PC scores omitted")

# ── Patient-level CAF2-AR+ proportion ──
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
meta_by_id = {f.stem: f for f in meta_files}
tumor_ids  = set(clinical.index[clinical.is_tumor == "yes"])
caf2_props: list[dict] = []
for sid in sorted(set(meta_by_id) & tumor_ids):
    mdf = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    if "sample_id" in mdf.index.names:
        mdf = mdf.xs(sid, level="sample_id") if sid in mdf.index.get_level_values("sample_id") else mdf
    lbl = mdf["meta_label"].values
    if len(lbl) == 0:
        continue
    caf2_props.append({"sample_id": sid, "prop": (lbl == "stromal-CAF2-AR+").sum() / len(lbl)})
caf2_df = pd.DataFrame(caf2_props).set_index("sample_id")
caf2_df["pat_id"] = pat_id_map(caf2_df.index)
scores["caf2_ar_plus_prop"] = caf2_df.groupby("pat_id")["prop"].mean()

# Normalize all series indices to string pat_ids from clinical
_float_to_str: dict[float, str] = {}
for p in clinical[clinical.is_tumor == "yes"]["pat_id"].unique():
    try:
        _float_to_str[float(p)] = p
    except ValueError:
        pass

def normalize_index(s: pd.Series) -> pd.Series:
    new_idx = [_float_to_str.get(x, x) if not isinstance(x, str) else x for x in s.index]
    return s.set_axis(new_idx)

scores = {k: normalize_index(v) for k, v in scores.items()}
print(f"  {len(scores)} scores collected: {list(scores.keys())}")

# %% ─── Assemble patient score matrix ─────────────────────────────────────────

score_df = pd.DataFrame(scores)
pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[
        ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time",
         "gleason_grp", "stromogenic_smc_loss_reactive_stroma_present", "inflammation"]
    ]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")
pat_data = score_df.join(pat_clin, how="inner")
score_cols = list(scores.keys())
print(f"  Patient matrix: {len(pat_data)} × {len(score_cols)} scores")

score_df.to_csv(OUT_TAB / "patient_all_scores.csv")

# %% ─── Spearman correlation matrix ──────────────────────────────────────────

print("Computing Spearman correlation matrix...")
corr_rows: list[dict] = []
for i, c1 in enumerate(score_cols):
    for c2 in score_cols[i:]:
        v1 = pat_data[c1].to_numpy(dtype=float) if isinstance(pat_data[c1], pd.Series) else pat_data[c1].iloc[:, 0].to_numpy(dtype=float)
        v2 = pat_data[c2].to_numpy(dtype=float) if isinstance(pat_data[c2], pd.Series) else pat_data[c2].iloc[:, 0].to_numpy(dtype=float)
        valid = ~(np.isnan(v1) | np.isnan(v2))
        if valid.sum() < 20:
            continue
        rho, p = spearmanr(v1[valid], v2[valid])
        corr_rows.append({"score_a": c1, "score_b": c2, "rho": float(rho), "p": float(p), "n": int(valid.sum())})

corr_pairs = pd.DataFrame(corr_rows)
corr_pairs.to_csv(OUT_TAB / "score_correlations.csv", index=False)

# Build symmetric matrix via explicit numpy assignment
_idx = {c: i for i, c in enumerate(score_cols)}
_mat = np.eye(len(score_cols))
for sa, sb, rho_v in zip(corr_pairs["score_a"], corr_pairs["score_b"], corr_pairs["rho"]):
    i, j = _idx[sa], _idx[sb]
    _mat[i, j] = float(rho_v)
    _mat[j, i] = float(rho_v)
rho_mat = pd.DataFrame(_mat, index=score_cols, columns=score_cols)

# %% ─── Figure 1: Correlation heatmap (clustered) ────────────────────────────

fig, ax = plt.subplots(figsize=(12, 10))
short_names = [s[:22] for s in score_cols]
g = sns.clustermap(
    rho_mat,
    cmap="RdBu_r", vmin=-1, vmax=1, center=0,
    xticklabels=short_names, yticklabels=short_names,
    figsize=(12, 10), annot=False,
    method="average",
)
g.ax_heatmap.set_title("Spearman correlation between patient-level TME scores")
g.fig.savefig(OUT_FIG / "score_correlation_clustermap.pdf", bbox_inches="tight")
plt.close()

# Non-clustered heatmap for clarity
fig, ax = plt.subplots(figsize=(11, 9))
mask = np.eye(len(score_cols), dtype=bool)
sns.heatmap(rho_mat, cmap="RdBu_r", vmin=-1, vmax=1, center=0,
            mask=mask, annot=True, fmt=".2f", ax=ax, annot_kws={"size": 6},
            xticklabels=short_names, yticklabels=short_names)
ax.set_title("Patient-level spatial score correlations (Spearman rho)")
plt.tight_layout()
fig.savefig(OUT_FIG / "score_correlation_heatmap.pdf", bbox_inches="tight")
plt.close()

# %% ─── Report strongest correlations ─────────────────────────────────────────

off_diag = corr_pairs[corr_pairs["score_a"] != corr_pairs["score_b"]].copy()
off_diag["abs_rho"] = off_diag["rho"].abs()
print("\nStrongest positive correlations:")
print(off_diag.nlargest(8, "rho")[["score_a", "score_b", "rho", "p", "n"]].to_string(index=False))
print("\nStrongest negative correlations:")
print(off_diag.nsmallest(8, "rho")[["score_a", "score_b", "rho", "p", "n"]].to_string(index=False))

# %% ─── Score correlations with clinical outcome ──────────────────────────────

print("\nCorrelating each score with Gleason and survival events...")
clin_corr_results: list[dict] = []
for col in score_cols:
    for clin_var in ["gleason_grp", "psa_progr", "clinical_progr"]:
        if clin_var not in pat_data.columns:
            continue
        v1 = pat_data[col].to_numpy(dtype=float, na_value=np.nan) if isinstance(pat_data[col], pd.Series) else pat_data[col].iloc[:, 0].to_numpy(dtype=float)
        v2 = pd.to_numeric(pat_data[clin_var], errors="coerce").to_numpy(dtype=float)
        valid = ~(np.isnan(v1) | np.isnan(v2))
        if valid.sum() < 20:
            continue
        rho, p = spearmanr(v1[valid], v2[valid])
        clin_corr_results.append({"score": col, "clinical": clin_var, "rho": float(rho), "p": float(p)})

clin_corr_df = pd.DataFrame(clin_corr_results)
_, fdr, _, _ = multipletests(clin_corr_df["p"].fillna(1), method="fdr_bh")
clin_corr_df["fdr"] = fdr
clin_corr_df.to_csv(OUT_TAB / "score_clinical_correlations.csv", index=False)

sig_clin = clin_corr_df[clin_corr_df["fdr"] < 0.1].sort_values("fdr")
print(f"  {len(sig_clin)} significant score-clinical associations (FDR<0.1):")
print(sig_clin[["score", "clinical", "rho", "p", "fdr"]].to_string(index=False))

# Heatmap: score × clinical variable (rho)
pivot_clin = clin_corr_df.pivot(index="score", columns="clinical", values="rho")
fig, ax = plt.subplots(figsize=(5, 10))
sns.heatmap(pivot_clin, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax,
            yticklabels=[s[:28] for s in pivot_clin.index])
ax.set_title("Spearman rho: TME scores vs clinical variables")
plt.tight_layout()
fig.savefig(OUT_FIG / "score_clinical_correlations.pdf", bbox_inches="tight")
plt.close()

# %% ─── Bivariate scatter: Shannon vs CAF1 periglandular (key pair) ───────────

for x_col, y_col, label in [
    ("h_global", "caf1_periglandular", "Shannon diversity vs CAF1 periglandular fraction"),
    ("h_global", "caf2_ar_plus_prop", "Shannon diversity vs CAF2-AR+ proportion"),
    ("caf1_periglandular", "caf2_ar_plus_prop", "CAF1 periglandular vs CAF2-AR+ proportion"),
]:
    if x_col not in pat_data.columns or y_col not in pat_data.columns:
        continue
    sub = pat_data[[x_col, y_col, "psa_progr", "gleason_grp"]].dropna()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, color_col, ctitle in zip(axes,
        ["gleason_grp", "psa_progr"],
        ["Gleason GG", "PSA recurrence"]):
        sc = ax.scatter(sub[x_col], sub[y_col], c=sub[color_col], cmap="RdYlBu_r",
                        s=20, alpha=0.7)
        plt.colorbar(sc, ax=ax, label=ctitle)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        rho, p = spearmanr(sub[x_col], sub[y_col])
        ax.set_title(f"rho={rho:.2f}, p={p:.3f}")
    fig.suptitle(label)
    plt.tight_layout()
    fname = f"scatter_{x_col}_vs_{y_col}.pdf".replace("/", "-")
    fig.savefig(OUT_FIG / fname, bbox_inches="tight")
    plt.close()

print("\nRQ15 complete. Outputs saved to", OUT_TAB, "and", OUT_FIG)
