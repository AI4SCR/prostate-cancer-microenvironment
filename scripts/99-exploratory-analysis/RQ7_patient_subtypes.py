"""
RQ7: Inter-patient Heterogeneity — Patient TME Subtypes

Clusters patients by mean ROI cell type composition, identifies stable patient subtypes,
and tests their associations with clinical variables and survival.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ7")
OUT_TAB = Path("output/tables/RQ7")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load and prepare data ─────────────────────────────────────────────────

print("Loading data...")
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
clinical = pd.read_parquet(BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet")
tumor_ids = set(clinical.index[clinical.is_tumor == "yes"])

cell_meta = pd.concat(
    [pd.read_parquet(f, engine="fastparquet") for f in meta_files], ignore_index=False
)
cell_meta = cell_meta.loc[cell_meta.index.get_level_values("sample_id").isin(tumor_ids)]

meta_counts = cell_meta.groupby(["sample_id", "meta_label"]).size().unstack(fill_value=0)
meta_props = meta_counts.div(meta_counts.sum(axis=1), axis=0)

# CLR transform
def clr(df: pd.DataFrame, pseudocount: float = 1e-5) -> pd.DataFrame:
    X = df.values + pseudocount
    log_X = np.log(X)
    return pd.DataFrame(log_X - log_X.mean(axis=1, keepdims=True), index=df.index, columns=df.columns)

props_clr = clr(meta_props)
props_clr_clin = props_clr.join(clinical[["pat_id"]])

ct_cols = list(meta_props.columns)

# %% ─── Patient-level mean CLR composition ────────────────────────────────────

pat_mean_clr = props_clr_clin.groupby("pat_id")[ct_cols].mean()
print(f"Patient-level profiles: {pat_mean_clr.shape}")

# Patient-level clinical
pat_clinical = clinical.groupby("pat_id").first()[
    ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time",
     "gleason_grp", "stromogenic_smc_loss_reactive_stroma_present",
     "inflammation", "cribriform", "d_amico_risk"]
]
pat_all = pat_mean_clr.join(pat_clinical)
pat_all.to_csv(OUT_TAB / "patient_mean_clr_composition.csv")

# %% ─── Variance decomposition ────────────────────────────────────────────────

print("Computing variance decomposition (between vs within patients)...")

# Total variance per cell type
total_var = props_clr.var()

# Between-patient variance: variance of patient means
between_var = pat_mean_clr.var()

# Within-patient variance: variance of ROI-level residuals (roi - patient mean)
roi_residuals = props_clr.join(clinical[["pat_id"]])
roi_residuals = roi_residuals[ct_cols] - roi_residuals.groupby("pat_id")[ct_cols].transform("mean")
within_var = roi_residuals.var()

var_decomp = pd.DataFrame({
    "total_var": total_var,
    "between_var": between_var,
    "within_var": within_var,
    "between_frac": between_var / (total_var + 1e-10),
    "within_frac": within_var / (total_var + 1e-10),
}).sort_values("between_frac", ascending=False)
var_decomp.to_csv(OUT_TAB / "variance_decomposition.csv")
print("Variance decomposition (sorted by between-patient fraction):")
print(var_decomp[["between_frac","within_frac"]].to_string())

# %% ─── Hierarchical clustering with stability assessment ─────────────────────

print("\nClustering patients by TME composition...")
np.random.seed(42)

# Standardize for clustering
X = StandardScaler().fit_transform(pat_mean_clr.values)
pat_ids = pat_mean_clr.index.tolist()

# Pairwise Euclidean distance
D = pdist(X, metric="euclidean")
Z = linkage(D, method="ward")

# Silhouette scores for k=2..6
print("Silhouette scores:")
sil_scores = {}
for k in range(2, 7):
    labels = fcluster(Z, k, criterion="maxclust")
    sil = silhouette_score(X, labels)
    sil_scores[k] = sil
    print(f"  k={k}: silhouette={sil:.4f}")

best_k = max(sil_scores, key=sil_scores.get)
print(f"\nBest k by silhouette: {best_k}")

# Bootstrap stability (100 iterations)
print("Assessing cluster stability (100 bootstraps)...")
stability_scores = {}
for k in [2, 3, 4]:
    agreements = []
    labels_full = fcluster(Z, k, criterion="maxclust")
    for _ in range(100):
        boot_idx = np.random.choice(len(X), len(X), replace=True)
        X_boot = X[boot_idx]
        D_boot = pdist(X_boot, metric="euclidean")
        Z_boot = linkage(D_boot, method="ward")
        labels_boot = fcluster(Z_boot, k, criterion="maxclust")
        # Match cluster labels via majority vote and compute agreement
        # Simplified: fraction of pairs with same relative cluster assignment
        n = len(boot_idx)
        same_orig = (labels_full[boot_idx][:, None] == labels_full[boot_idx][None, :])
        same_boot = (labels_boot[:, None] == labels_boot[None, :])
        agreement = (same_orig == same_boot).mean()
        agreements.append(agreement)
    stability_scores[k] = np.mean(agreements)
    print(f"  k={k}: mean pairwise agreement={np.mean(agreements):.4f}")

# Use k=3 as a biologically interpretable choice
K = 3
cluster_labels = fcluster(Z, K, criterion="maxclust")
pat_all["cluster"] = cluster_labels
pat_all["cluster"] = pat_all["cluster"].astype(str)
pat_all.to_csv(OUT_TAB / "patient_clusters.csv")
print(f"\nCluster sizes (k={K}):")
print(pat_all["cluster"].value_counts().sort_index().to_dict())

# %% ─── Figure 1: Dendrogram + heatmap ────────────────────────────────────────

fig = plt.figure(figsize=(18, 10))
ax_dendro = fig.add_axes([0.05, 0.1, 0.15, 0.8])
ax_heat   = fig.add_axes([0.22, 0.1, 0.65, 0.8])
ax_cbar   = fig.add_axes([0.89, 0.1, 0.02, 0.8])

# Dendrogram
dend = dendrogram(Z, ax=ax_dendro, orientation="left",
                  no_labels=True, color_threshold=Z[-K+1, 2])
ax_dendro.set_xlabel("Distance")
ax_dendro.set_title("Ward linkage")

# Reorder heatmap by dendrogram
order = dend["leaves"][::-1]
X_ordered = pat_mean_clr.values[order]
cluster_ordered = cluster_labels[order]

im = ax_heat.imshow(X_ordered, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax_heat.set_xticks(range(len(ct_cols)))
ax_heat.set_xticklabels(ct_cols, rotation=45, ha="right", fontsize=8)
ax_heat.set_yticks([])
ax_heat.set_title("Patient TME composition (CLR, z-scored)\ncolumns = cell types, rows = patients")

# Cluster color bar on right
cmap_clusters = plt.cm.Set1
for i, cl in enumerate(cluster_ordered):
    ax_heat.add_patch(plt.Rectangle(
        (len(ct_cols) - 0.5, i - 0.5), 0.8, 1,
        color=cmap_clusters((cl - 1) / K), clip_on=False
    ))

plt.colorbar(im, cax=ax_cbar, label="CLR proportion (z-scored)")
fig.savefig(OUT_FIG / "patient_clustering_heatmap.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 2: PCA colored by cluster and clinical variables ───────────────

pca = PCA(n_components=2)
pca_coords = pca.fit_transform(X)
pca_df = pd.DataFrame(pca_coords, index=pat_ids, columns=["PC1","PC2"])
pca_df = pca_df.join(pat_all[["cluster","gleason_grp","psa_progr","clinical_progr",
                                "stromogenic_smc_loss_reactive_stroma_present"]])

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
palette = {"1": "#e41a1c", "2": "#377eb8", "3": "#4daf4a", "4": "#984ea3", "5": "#ff7f00"}

ax = axes[0]
for cl, sub in pca_df.groupby("cluster"):
    ax.scatter(sub.PC1, sub.PC2, c=palette.get(str(cl), "gray"), s=30, alpha=0.7, label=f"Cluster {cl}")
ax.set_title(f"Patient clusters (k={K})")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.legend(fontsize=8)

ax = axes[1]
sc = ax.scatter(pca_df.PC1, pca_df.PC2, c=pca_df.gleason_grp, cmap="RdYlBu_r", s=30, alpha=0.7)
plt.colorbar(sc, ax=ax, label="Gleason grade group")
ax.set_title("Colored by Gleason grade group")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")

ax = axes[2]
colors_psa = pca_df["psa_progr"].map({0: "steelblue", 1: "red"}).fillna("gray")
ax.scatter(pca_df.PC1, pca_df.PC2, c=colors_psa, s=30, alpha=0.7)
ax.set_title("Colored by PSA recurrence\n(blue=no, red=yes)")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")

plt.tight_layout()
fig.savefig(OUT_FIG / "patient_pca_clusters.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 3: Cluster profiles ───────────────────────────────────────────

cluster_means = pat_all.groupby("cluster")[ct_cols].mean()
cluster_means_z = (cluster_means - cluster_means.mean()) / (cluster_means.std() + 1e-10)

fig, ax = plt.subplots(figsize=(13, 5))
sns.heatmap(cluster_means_z.T, cmap="RdBu_r", center=0, linewidths=0.5, ax=ax,
            cbar_kws={"label": "z-score (CLR proportion)"},
            annot=True, fmt=".2f", annot_kws={"size": 7})
ax.set_title("Patient cluster profiles — mean CLR composition (z-scored across clusters)")
ax.set_xlabel("Cluster")
ax.set_ylabel("Cell type (meta-label)")
plt.tight_layout()
fig.savefig(OUT_FIG / "cluster_profiles_heatmap.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 4: Cluster vs clinical variables ───────────────────────────────

CLIN_VARS = ["gleason_grp", "stromogenic_smc_loss_reactive_stroma_present",
              "inflammation", "cribriform", "psa_progr", "clinical_progr"]
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()

for ax, col in zip(axes, CLIN_VARS):
    sub = pat_all.dropna(subset=["cluster", col])
    if sub[col].dtype == object or sub[col].nunique() <= 5:
        ct_tab = pd.crosstab(sub["cluster"], sub[col], normalize="index") * 100
        ct_tab.plot(kind="bar", ax=ax, stacked=True, colormap="Set2", legend=True)
        ax.set_title(f"{col} distribution per cluster")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("% patients")
        ax.tick_params(axis="x", rotation=0)
    else:
        sns.boxplot(data=sub, x="cluster", y=col, ax=ax, palette="Set1")
        ax.set_title(f"{col} by cluster")

plt.tight_layout()
fig.savefig(OUT_FIG / "cluster_clinical_distributions.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Cluster-specific Kruskal-Wallis for each cell type ───────────────────

print("\nKruskal-Wallis test per cell type across clusters:")
kw_results = []
for ct in ct_cols:
    groups = [pat_all[pat_all.cluster == cl][ct].dropna().values
              for cl in sorted(pat_all.cluster.unique())]
    groups = [g for g in groups if len(g) >= 3]
    if len(groups) < 2:
        continue
    try:
        stat, pval = kruskal(*groups)
        kw_results.append({"cell_type": ct, "stat": stat, "pval": pval})
    except ValueError:
        pass

kw_df = pd.DataFrame(kw_results)
kw_df["fdr"] = multipletests(kw_df["pval"].fillna(1), method="fdr_bh")[1]
kw_df = kw_df.sort_values("fdr")
kw_df.to_csv(OUT_TAB / "cluster_kruskal_wallis.csv", index=False)
print(f"Significant (FDR < 0.1): {(kw_df.fdr < 0.1).sum()}")
print(kw_df[kw_df.fdr < 0.1][["cell_type","fdr"]].to_string(index=False))

# %% ─── Survival analysis by cluster ─────────────────────────────────────────

print("\nKaplan-Meier by cluster...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
cluster_colors = {"1": "#e41a1c", "2": "#377eb8", "3": "#4daf4a"}

for ax, (event_col, time_col, label) in zip(axes, [
    ("psa_progr", "psa_progr_time", "PSA recurrence"),
    ("clinical_progr", "clinical_progr_time", "Clinical progression"),
]):
    sub = pat_all.dropna(subset=["cluster", event_col, time_col])
    clusters = sorted(sub["cluster"].unique())

    for cl in clusters:
        g = sub[sub.cluster == cl]
        KaplanMeierFitter().fit(
            g[time_col], g[event_col], label=f"Cluster {cl} (n={len(g)})"
        ).plot_survival_function(ax=ax, ci_show=False, color=cluster_colors.get(cl, "gray"))

    # Log-rank across all clusters
    groups_data = [(sub[sub.cluster == cl][time_col], sub[sub.cluster == cl][event_col])
                   for cl in clusters]
    if len(groups_data) >= 2:
        from lifelines.statistics import multivariate_logrank_test
        mlr = multivariate_logrank_test(
            sub[time_col], sub["cluster"], sub[event_col]
        )
        ax.set_title(f"KM by cluster — {label}\nmultivariate log-rank p={mlr.p_value:.3f}")
    else:
        ax.set_title(f"KM by cluster — {label}")

    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.legend(fontsize=8)

plt.tight_layout()
fig.savefig(OUT_FIG / "km_by_cluster.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure: Variance decomposition bar chart ─────────────────────────────

fig, ax = plt.subplots(figsize=(10, 6))
var_decomp_sorted = var_decomp.sort_values("between_frac", ascending=True)
x = range(len(var_decomp_sorted))
ax.barh(var_decomp_sorted.index, var_decomp_sorted["between_frac"],
        color="steelblue", alpha=0.8, label="Between-patient")
ax.barh(var_decomp_sorted.index, var_decomp_sorted["within_frac"],
        left=var_decomp_sorted["between_frac"],
        color="salmon", alpha=0.8, label="Within-patient")
ax.axvline(0.5, color="k", linestyle="--", linewidth=0.8, alpha=0.5)
ax.set_xlabel("Fraction of total variance")
ax.set_title("Variance decomposition: between- vs within-patient\n(CLR-transformed cell type proportions)")
ax.legend()
plt.tight_layout()
fig.savefig(OUT_FIG / "variance_decomposition.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Summary ───────────────────────────────────────────────────────────────

print(f"\n=== RQ7 Summary ===")
print(f"Patients clustered: {len(pat_all)}")
print(f"\nCluster sizes: {pat_all.cluster.value_counts().sort_index().to_dict()}")
print(f"\nSilhouette scores: {sil_scores}")
print(f"Bootstrap stability: {stability_scores}")
print(f"\nVariance decomposition (between-patient fraction, top 5):")
print(var_decomp["between_frac"].head(5).to_string())
print(f"\nOutputs → {OUT_TAB}/ and {OUT_FIG}/")
