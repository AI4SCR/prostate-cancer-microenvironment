"""
RQ19: Cribriform Spatial Signature

Identifies spatial features that distinguish cribriform from non-cribriform ROIs
beyond cell-type composition. Uses:
  1. Epithelial self-contact enrichment (luminal clustering) — RQ6 showed this drops in cribriform
  2. CAF2-AR+ depletion (strongest compositional cribriform signal from RQ1)
  3. Spatial regularity of luminal epithelium (NN distance SD as a proxy for gland disorganization)
  4. Multi-feature classifier: logistic regression on spatial features predicting cribriform label

Clinically relevant because:
  - Cribriform is adversely prognostic and now guideline-relevant (equivalent to GG4)
  - Pathologist detection is subjective and inter-observer variable
  - A computational spatial signature could support automated cribriform detection
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, RocCurveDisplay
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ19")
OUT_TAB = Path("output/tables/RQ19")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32
N_PERM = 20

LUMINAL_LABELS = {"epithelial-luminal", "epithelial-luminal(Ki67+)"}
TUMOR_LABELS   = {"epithelial-(ERG+CD44+)", "epithelial-luminal(ERG+)",
                  "epithelial-luminal(ERG+p53+)"}
ALL_EPI_LABELS = LUMINAL_LABELS | TUMOR_LABELS

CAF2_AR_LABEL  = "stromal-CAF2(AR+)"
CAF2_LABELS    = {"stromal-CAF2(AR+)", "stromal-CAF2(AR+CES1+)", "stromal-CAF2(AR+EGR1+)"}
CAF1_PLUS      = "stromal-CAF1(CD105+)"

# %% ─── Load data ─────────────────────────────────────────────────────────────

print("Loading data...")
meta_files    = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
spatial_files = sorted((BASE_DIR / "02_processed/features/spatial/filtered-annotated").glob("*.parquet"))
meta_by_id    = {f.stem: f for f in meta_files}
spatial_by_id = {f.stem: f for f in spatial_files}

clinical = pd.read_parquet(
    BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet"
)
tumor_ids  = set(clinical.index[clinical.is_tumor == "yes"])
sample_ids = sorted(set(meta_by_id) & set(spatial_by_id) & tumor_ids)

# Keep only ROIs with known cribriform annotation
crib_known = clinical.loc[clinical.index.isin(sample_ids), "cribriform"].isin(["yes", "no"])
sample_ids = [sid for sid in sample_ids if crib_known.get(sid, False)]
print(f"  {len(sample_ids)} tumor ROIs with cribriform annotation")
n_crib = clinical.loc[sample_ids, "cribriform"].eq("yes").sum()
print(f"  Cribriform: {n_crib}, Non-cribriform: {len(sample_ids) - n_crib}")


def contact_enrichment(coords: np.ndarray, mask_a: np.ndarray, mask_b: np.ndarray,
                       radius: float, n_perm: int = 20) -> float:
    """Log2(observed/expected) contact enrichment for cell type A vs B."""
    tree = cKDTree(coords)
    a_idx = np.where(mask_a)[0]
    b_idx = np.where(mask_b)[0]
    if len(a_idx) == 0 or len(b_idx) == 0:
        return 0.0
    # Observed contacts
    obs = sum(
        1 for ai in a_idx
        for ni in tree.query_ball_point(coords[ai], r=radius)
        if ni != ai and mask_b[ni]
    )
    # Expected via permutation
    labels = np.zeros(len(coords), dtype=int)
    labels[mask_a] = 1
    labels[mask_b] = 2
    exp_counts = []
    rng = np.random.default_rng(42)
    for _ in range(n_perm):
        perm = labels.copy()
        rng.shuffle(perm)
        a_p = np.where(perm == 1)[0]
        b_p = set(np.where(perm == 2)[0])
        exp = sum(
            1 for ai in a_p
            for ni in tree.query_ball_point(coords[ai], r=radius)
            if ni != ai and ni in b_p
        )
        exp_counts.append(exp)
    expected = np.mean(exp_counts)
    return float(np.log2((obs + 0.5) / (expected + 0.5)))


# %% ─── Per-ROI spatial features for cribriform signature ─────────────────────

print("Computing spatial features per ROI...")
roi_records: list[dict] = []

for sid in sample_ids:
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    spat_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")

    if "sample_id" in meta_df.index.names:
        meta_df = meta_df.xs(sid, level="sample_id") if sid in meta_df.index.get_level_values("sample_id") else meta_df
    if "sample_id" in spat_df.index.names:
        spat_df = spat_df.xs(sid, level="sample_id") if sid in spat_df.index.get_level_values("sample_id") else spat_df

    idx = meta_df.index.intersection(spat_df.index)
    labels = meta_df.loc[idx, "label"].to_numpy()
    coords = spat_df.loc[idx, ["x", "y"]].to_numpy(dtype=np.float32)
    n = len(labels)

    epi_mask  = np.isin(labels, list(ALL_EPI_LABELS))
    lum_mask  = np.isin(labels, list(LUMINAL_LABELS))
    caf2_mask = np.isin(labels, list(CAF2_LABELS))
    caf1_mask = labels == CAF1_PLUS

    # Cell type proportions
    caf2_prop = caf2_mask.sum() / n if n > 0 else 0.0
    caf1_prop = caf1_mask.sum() / n if n > 0 else 0.0
    epi_prop  = epi_mask.sum() / n if n > 0 else 0.0

    # Luminal epithelium spatial regularity:
    # NN distance among luminal cells — high SD = irregular gland arrangement
    lum_coords = coords[lum_mask]
    lum_nn_dist_mean = np.nan
    lum_nn_dist_sd   = np.nan
    if len(lum_coords) >= 5:
        tree_lum = cKDTree(lum_coords)
        dists, _ = tree_lum.query(lum_coords, k=2)  # k=2: nearest neighbor (not self)
        nn_dists = dists[:, 1]
        lum_nn_dist_mean = float(nn_dists.mean())
        lum_nn_dist_sd   = float(nn_dists.std())

    # Luminal self-contact enrichment (log2 obs/expected)
    # proxy for how tightly luminal cells cluster
    lum_self_enrichment = np.nan
    if lum_mask.sum() >= 5:
        lum_self_enrichment = contact_enrichment(coords, lum_mask, lum_mask, RADIUS, N_PERM)

    # CAF2-AR+ exclusion from luminal epithelium
    caf2_lum_enrichment = np.nan
    if caf2_mask.sum() >= 3 and lum_mask.sum() >= 3:
        caf2_lum_enrichment = contact_enrichment(coords, caf2_mask, lum_mask, RADIUS, N_PERM)

    # CAF1-CD105+ contact with luminal epithelium
    caf1_lum_enrichment = np.nan
    if caf1_mask.sum() >= 3 and lum_mask.sum() >= 3:
        caf1_lum_enrichment = contact_enrichment(coords, caf1_mask, lum_mask, RADIUS, N_PERM)

    roi_records.append({
        "sample_id": sid,
        "n": n,
        "caf2_prop": caf2_prop,
        "caf1_prop": caf1_prop,
        "epi_prop": epi_prop,
        "lum_nn_dist_mean": lum_nn_dist_mean,
        "lum_nn_dist_sd": lum_nn_dist_sd,
        "lum_self_enrichment": lum_self_enrichment,
        "caf2_lum_enrichment": caf2_lum_enrichment,
        "caf1_lum_enrichment": caf1_lum_enrichment,
    })

roi_df = pd.DataFrame(roi_records).set_index("sample_id")
roi_df["cribriform"] = clinical.loc[roi_df.index, "cribriform"]
roi_df["gleason_grp"] = pd.to_numeric(clinical.loc[roi_df.index, "gleason_grp"], errors="coerce")
roi_df.to_csv(OUT_TAB / "roi_cribriform_spatial_features.csv")
print(f"  {len(roi_df)} ROIs computed")

# %% ─── Univariate associations with cribriform ───────────────────────────────

print("\nUnivariate associations with cribriform label...")
spatial_cols = ["caf2_prop", "caf1_prop", "epi_prop",
                "lum_nn_dist_mean", "lum_nn_dist_sd",
                "lum_self_enrichment", "caf2_lum_enrichment", "caf1_lum_enrichment"]

results: list[dict] = []
for col in spatial_cols:
    crib_vals = roi_df.loc[roi_df["cribriform"] == "yes", col].dropna()
    non_vals  = roi_df.loc[roi_df["cribriform"] == "no",  col].dropna()
    if len(crib_vals) < 5 or len(non_vals) < 5:
        continue
    stat, p = mannwhitneyu(crib_vals, non_vals, alternative="two-sided")
    rbc = 1 - 2 * stat / (len(crib_vals) * len(non_vals))
    results.append({"feature": col, "rbc": rbc, "p": p,
                    "n_crib": len(crib_vals), "n_non": len(non_vals),
                    "mean_crib": crib_vals.mean(), "mean_non": non_vals.mean()})

res_df = pd.DataFrame(results)
_, fdr, _, _ = multipletests(res_df["p"].fillna(1), method="fdr_bh")
res_df["fdr"] = fdr
res_df = res_df.sort_values("p")
res_df.to_csv(OUT_TAB / "cribriform_univariate.csv", index=False)
print(res_df[["feature","rbc","p","fdr","mean_crib","mean_non"]].to_string(index=False))

# %% ─── Multi-feature logistic regression: cribriform classifier ──────────────

print("\nMulti-feature logistic regression for cribriform detection...")
feature_cols = [c for c in spatial_cols if roi_df[c].notna().sum() > 100]
sub = roi_df[feature_cols + ["cribriform"]].dropna()
y = (sub["cribriform"] == "yes").astype(int)
X_raw = sub[feature_cols].values

scaler = StandardScaler()
X = scaler.fit_transform(X_raw)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
lr = LogisticRegression(C=1.0, max_iter=500)
aucs = cross_val_score(lr, X, y, cv=cv, scoring="roc_auc")
print(f"  Features: {feature_cols}")
print(f"  5-fold CV AUC: {aucs.mean():.3f} ± {aucs.std():.3f}")

# Full model for coefficients
lr.fit(X, y)
coef_df = pd.DataFrame({
    "feature": feature_cols,
    "coef": lr.coef_[0],
    "abs_coef": np.abs(lr.coef_[0]),
}).sort_values("abs_coef", ascending=False)
coef_df.to_csv(OUT_TAB / "cribriform_classifier_coefs.csv", index=False)
print(coef_df[["feature", "coef"]].to_string(index=False))

# Compare: composition-only AUC vs spatial+composition AUC
comp_only_cols = [c for c in ["caf2_prop", "caf1_prop", "epi_prop"] if c in feature_cols]
X_comp = scaler.fit_transform(sub[comp_only_cols].values)
aucs_comp = cross_val_score(LogisticRegression(C=1.0, max_iter=500), X_comp, y,
                             cv=cv, scoring="roc_auc")
print(f"\n  Composition-only AUC: {aucs_comp.mean():.3f} ± {aucs_comp.std():.3f}")
print(f"  Spatial+composition AUC: {aucs.mean():.3f} ± {aucs.std():.3f}")
print(f"  Delta AUC: +{aucs.mean() - aucs_comp.mean():.3f}")

# %% ─── Figures ───────────────────────────────────────────────────────────────

# 1. Violin plots: spatial features by cribriform
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
for ax, col in zip(axes.ravel(), spatial_cols):
    sub_plot = roi_df[["cribriform", col]].dropna()
    sns.violinplot(data=sub_plot, x="cribriform", y=col, order=["no", "yes"],
                   ax=ax, inner="box", cut=0)
    ax.set_title(col[:30], fontsize=8)
    ax.set_xlabel("Cribriform")
    row = res_df[res_df["feature"] == col]
    if not row.empty:
        ax.set_title(f"{col[:26]}\nrbc={row['rbc'].values[0]:.2f}, FDR={row['fdr'].values[0]:.3f}", fontsize=7)
plt.suptitle("Spatial features: cribriform vs non-cribriform")
plt.tight_layout()
fig.savefig(OUT_FIG / "spatial_features_cribriform.pdf", bbox_inches="tight")
plt.close()

# 2. ROC curve: full model vs composition-only
from sklearn.model_selection import cross_val_predict
y_score_full = cross_val_predict(LogisticRegression(C=1.0, max_iter=500),
                                  X, y, cv=cv, method="predict_proba")[:, 1]
y_score_comp = cross_val_predict(LogisticRegression(C=1.0, max_iter=500),
                                  X_comp, y, cv=cv, method="predict_proba")[:, 1]

from sklearn.metrics import roc_curve
fig, ax = plt.subplots(figsize=(5, 5))
for y_sc, label, color in [
    (y_score_full, f"Spatial+composition (AUC={roc_auc_score(y, y_score_full):.3f})", "steelblue"),
    (y_score_comp, f"Composition only (AUC={roc_auc_score(y, y_score_comp):.3f})", "tomato"),
]:
    fpr, tpr, _ = roc_curve(y, y_sc)
    ax.plot(fpr, tpr, label=label, color=color)
ax.plot([0, 1], [0, 1], "k--", linewidth=0.8)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("Cribriform detection: spatial vs composition features\n5-fold CV ROC")
ax.legend(fontsize=8)
plt.tight_layout()
fig.savefig(OUT_FIG / "roc_cribriform_classifier.pdf", bbox_inches="tight")
plt.close()

print("\nRQ19 complete. Outputs →", OUT_TAB, "and", OUT_FIG)
