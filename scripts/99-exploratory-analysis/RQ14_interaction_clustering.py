"""
RQ14: Patient-Level Clustering on Spatial Interaction Profiles

Tests whether patients cluster more distinctly based on their spatial wiring
(17×17 cell-cell interaction enrichment matrix from RQ6) than on composition alone (RQ7).
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ14")
OUT_TAB = Path("output/tables/RQ14")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load per-ROI interaction enrichment ───────────────────────────────────

print("Loading per-ROI interaction enrichment matrix...")
roi_int = pd.read_parquet("output/tables/RQ6/roi_interaction_enrichment.parquet")
print(f"  {roi_int.shape[0]} ROIs × {roi_int.shape[1]} interaction features")

# %% ─── Load clinical data ────────────────────────────────────────────────────

clinical = pd.read_parquet(
    BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet"
)
pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[
        ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time",
         "gleason_grp", "stromogenic_smc_loss_reactive_stroma_present", "inflammation"]
    ]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")

# %% ─── Aggregate to patient level ────────────────────────────────────────────

roi_int.index = roi_int.index.astype(str)
roi_int["pat_id"] = clinical.loc[roi_int.index.intersection(clinical.index), "pat_id"]
int_cols = [c for c in roi_int.columns if c != "pat_id"]

pat_int = roi_int.groupby("pat_id")[int_cols].mean()
print(f"  {len(pat_int)} patients with interaction profiles")

# Impute missing (NaN from ROIs where cell types absent) with 0
pat_int = pat_int.fillna(0.0)

# %% ─── Hierarchical clustering ───────────────────────────────────────────────

print("Hierarchical clustering on patient interaction profiles...")
X = pat_int.values

# Silhouette scores for k=2..6
sil_scores: dict[int, float] = {}
Z = linkage(X, method="ward", metric="euclidean")
for k in range(2, 7):
    labels = fcluster(Z, k, criterion="maxclust")
    if len(np.unique(labels)) > 1:
        sil_scores[k] = silhouette_score(X, labels)

print(f"  Silhouette scores: {sil_scores}")
best_k = max(sil_scores, key=sil_scores.get)
print(f"  Best k={best_k} (silhouette={sil_scores[best_k]:.3f})")

# Compare with RQ7 composition-based silhouette
print("  [RQ7 reference: max silhouette ≈ 0.12 at k=2 on composition]")

# Assign clusters at k=2 and k=best_k
for k in [2, best_k]:
    cluster_labels = fcluster(Z, k, criterion="maxclust")
    pat_int[f"cluster_k{k}"] = cluster_labels

pat_int.to_csv(OUT_TAB / "patient_interaction_clusters.csv")

sil_df = pd.Series(sil_scores, name="silhouette").reset_index().rename(columns={"index": "k"})
sil_df.to_csv(OUT_TAB / "silhouette_scores.csv", index=False)

# %% ─── Figure 1: Silhouette scores ──────────────────────────────────────────

fig, ax = plt.subplots(figsize=(5, 3))
ax.plot(list(sil_scores.keys()), list(sil_scores.values()), "o-")
ax.set_xlabel("k (number of clusters)")
ax.set_ylabel("Silhouette score")
ax.set_title("Patient clustering on interaction profiles\nvs k=2 benchmark (RQ7 composition: 0.12)")
ax.axhline(0.12, color="gray", linestyle="--", label="RQ7 composition benchmark")
ax.legend(fontsize=8)
plt.tight_layout()
fig.savefig(OUT_FIG / "silhouette_scores.pdf", bbox_inches="tight")
plt.close()

# %% ─── Figure 2: Dendrogram ──────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(12, 4))
dendrogram(Z, ax=ax, no_labels=True, color_threshold=0.7 * max(Z[:, 2]))
ax.set_title("Patient dendrogram — interaction enrichment profiles")
plt.tight_layout()
fig.savefig(OUT_FIG / "dendrogram.pdf", bbox_inches="tight")
plt.close()

# %% ─── Figure 3: PCA of interaction space colored by cluster ────────────────

pca = PCA(n_components=2)
coords = pca.fit_transform(X)
pca_df = pd.DataFrame(coords, index=pat_int.index, columns=["PC1", "PC2"])
pca_df[f"cluster_k{best_k}"] = pat_int[f"cluster_k{best_k}"]

fig, ax = plt.subplots(figsize=(6, 5))
for cl in sorted(pca_df[f"cluster_k{best_k}"].unique()):
    m = pca_df[pca_df[f"cluster_k{best_k}"] == cl]
    ax.scatter(m["PC1"], m["PC2"], label=f"C{cl}", s=20, alpha=0.8)
ax.legend()
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
ax.set_title(f"Patient interaction profiles — k={best_k} clusters")
plt.tight_layout()
fig.savefig(OUT_FIG / f"pca_interaction_k{best_k}.pdf", bbox_inches="tight")
plt.close()

# %% ─── Clinical associations ─────────────────────────────────────────────────

print("\nClinical associations with interaction-based clusters...")
pat_data = pat_int[[f"cluster_k{best_k}"]].join(pat_clin, how="inner")
print(f"  {len(pat_data)} patients with clinical data")

# Cluster sizes and clinical variable distributions
for var in ["gleason_grp", "stromogenic_smc_loss_reactive_stroma_present", "inflammation"]:
    if var not in pat_data.columns:
        continue
    ct = pd.crosstab(pat_data[f"cluster_k{best_k}"], pat_data[var])
    print(f"\n  {var}:")
    print(ct)

# %% ─── Cox PH: cluster k=best_k vs survival ─────────────────────────────────

print("\nCox PH: interaction clusters vs survival...")
cox_results: list[dict] = []
for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    sub = pat_data[[duration_col, event_col, f"cluster_k{best_k}"]].dropna()
    sub = sub[sub[duration_col] > 0]
    # Use cluster as categorical dummy — compare each cluster to cluster 1
    sub_dummy = pd.get_dummies(sub[f"cluster_k{best_k}"].astype(str), prefix="cl", drop_first=True)
    cph_df = pd.concat([sub[[duration_col, event_col]], sub_dummy], axis=1)
    cph_df = cph_df.rename(columns={duration_col: "T", event_col: "E"})
    if cph_df.shape[1] <= 2:
        continue
    try:
        cph = CoxPHFitter()
        cph.fit(cph_df, duration_col="T", event_col="E", show_progress=False)
        for coef in cph.summary.index:
            s = cph.summary.loc[coef]
            cox_results.append({
                "cluster_coef": coef, "endpoint": endpoint,
                "hr": np.exp(s["coef"]),
                "ci_low": np.exp(s["coef lower 95%"]),
                "ci_high": np.exp(s["coef upper 95%"]),
                "p": s["p"],
            })
    except Exception as e:
        print(f"  Cox failed: {e}")

cox_df = pd.DataFrame(cox_results)
if len(cox_df):
    cox_df.to_csv(OUT_TAB / "cluster_cox_results.csv", index=False)
    print(cox_df[["cluster_coef", "endpoint", "hr", "ci_low", "ci_high", "p"]].to_string(index=False))

# %% ─── KM curves: best k ─────────────────────────────────────────────────────

for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
]:
    sub = pat_data[[duration_col, event_col, f"cluster_k{best_k}"]].dropna()
    sub = sub[sub[duration_col] > 0]
    clusters = sorted(sub[f"cluster_k{best_k}"].unique())

    fig, ax = plt.subplots(figsize=(6, 4))
    for cl in clusters:
        m = sub[sub[f"cluster_k{best_k}"] == cl]
        KaplanMeierFitter().fit(m[duration_col], m[event_col], label=f"C{cl} (n={len(m)})").plot_survival_function(ax=ax)

    if len(clusters) == 2:
        c1, c2 = clusters
        m1, m2 = sub[sub[f"cluster_k{best_k}"] == c1], sub[sub[f"cluster_k{best_k}"] == c2]
        lr = logrank_test(m1[duration_col], m2[duration_col], m1[event_col], m2[event_col])
        ax.set_title(f"KM — {endpoint}\nInteraction clusters k={best_k} (log-rank p={lr.p_value:.3f})")
    else:
        ax.set_title(f"KM — {endpoint}\nInteraction clusters k={best_k}")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"km_{endpoint}_k{best_k}.pdf", bbox_inches="tight")
    plt.close()

print("\nRQ14 complete. Outputs saved to", OUT_TAB, "and", OUT_FIG)
