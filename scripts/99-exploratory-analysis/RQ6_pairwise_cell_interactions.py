"""
RQ6: Full Pairwise Cell-Cell Interaction Profiles

Computes permutation-based enrichment scores for all meta-label × meta-label contact pairs
in radius-32 neighborhoods, identifies the average PCa TME wiring diagram, and tests
which interactions differ across Gleason grade and other clinical groups.
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
from scipy.stats import mannwhitneyu, spearmanr, kruskal
from statsmodels.stats.multitest import multipletests
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ6")
OUT_TAB = Path("output/tables/RQ6")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32
N_PERM = 30  # permutations per ROI for enrichment

# %% ─── Load data ─────────────────────────────────────────────────────────────

print("Loading data...")
meta_files   = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
spatial_files = sorted((BASE_DIR / "02_processed/features/spatial/filtered-annotated").glob("*.parquet"))
meta_by_id   = {f.stem: f for f in meta_files}
spatial_by_id = {f.stem: f for f in spatial_files}

clinical = pd.read_parquet(BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet")
tumor_ids = set(clinical.index[clinical.is_tumor == "yes"])
sample_ids = sorted(set(meta_by_id) & set(spatial_by_id) & tumor_ids)
print(f"  {len(sample_ids)} tumor ROIs")

# Get full meta_label set
all_labels: set = set()
for sid in sample_ids[:30]:
    m = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    m = m.loc[sid] if sid in m.index.get_level_values("sample_id") else m
    all_labels.update(m["meta_label"].unique())
LABELS = sorted(all_labels)
L = len(LABELS)
label_idx = {l: i for i, l in enumerate(LABELS)}
print(f"  {L} meta-labels")

# %% ─── Per-ROI interaction enrichment ────────────────────────────────────────

print("Computing pairwise contact enrichment per ROI...")


def interaction_enrichment(coords: np.ndarray, labels: np.ndarray,
                             label_idx: dict, n_labels: int,
                             radius: float, n_perm: int) -> np.ndarray:
    """
    Returns log2(obs/exp + eps) enrichment matrix [n_labels × n_labels].
    obs[i,j] = number of radius contacts from cell type i to cell type j.
    exp = mean over permutations of shuffled labels.
    """
    n = len(labels)
    tree = cKDTree(coords)

    # Observed contact matrix
    obs = np.zeros((n_labels, n_labels), dtype=float)
    for i in range(n):
        neighbors = tree.query_ball_point(coords[i], r=radius)
        li = label_idx.get(labels[i], -1)
        if li < 0:
            continue
        for j in neighbors:
            if j == i:
                continue
            lj = label_idx.get(labels[j], -1)
            if lj >= 0:
                obs[li, lj] += 1

    # Expected: permuted label matrices
    exp_sum = np.zeros((n_labels, n_labels), dtype=float)
    perm_labels = labels.copy()
    for _ in range(n_perm):
        np.random.shuffle(perm_labels)
        perm_mat = np.zeros((n_labels, n_labels), dtype=float)
        for i in range(n):
            neighbors = tree.query_ball_point(coords[i], r=radius)
            li = label_idx.get(perm_labels[i], -1)
            if li < 0:
                continue
            for j in neighbors:
                if j == i:
                    continue
                lj = label_idx.get(perm_labels[j], -1)
                if lj >= 0:
                    perm_mat[li, lj] += 1
        exp_sum += perm_mat
    exp = exp_sum / n_perm

    # log2 enrichment
    eps = 0.5
    enrich = np.log2((obs + eps) / (exp + eps))
    return enrich


np.random.seed(42)
roi_enrichments = []

for i, sid in enumerate(sample_ids):
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(sample_ids)}")

    meta_df    = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    spatial_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")
    meta_df    = meta_df.loc[sid] if sid in meta_df.index.get_level_values("sample_id") else meta_df
    spatial_df = spatial_df.loc[sid] if sid in spatial_df.index.get_level_values("sample_id") else spatial_df
    common = meta_df.index.intersection(spatial_df.index)
    if len(common) < 20:
        continue
    meta_df    = meta_df.loc[common]
    spatial_df = spatial_df.loc[common]

    labels = meta_df["meta_label"].values
    coords = spatial_df[["x", "y"]].values

    enrich = interaction_enrichment(coords, labels, label_idx, L, RADIUS, N_PERM)

    # Flatten to named columns
    row = {"sample_id": sid}
    for a, la in enumerate(LABELS):
        for b, lb in enumerate(LABELS):
            row[f"{la}__x__{lb}"] = enrich[a, b]
    roi_enrichments.append(row)

enrich_df = pd.DataFrame(roi_enrichments).set_index("sample_id")
print(f"  {len(enrich_df)} ROIs processed; {enrich_df.shape[1]} interaction pairs")
enrich_df.to_parquet(OUT_TAB / "roi_interaction_enrichment.parquet")

# Merge with clinical
enrich_clin = enrich_df.join(clinical[["pat_id", "gleason_grp",
                                         "stromogenic_smc_loss_reactive_stroma_present",
                                         "inflammation", "cribriform",
                                         "psa_progr", "clinical_progr"]])

# %% ─── Mean interaction matrix (average PCa TME wiring) ────────────────────

print("Computing mean interaction matrix...")
interaction_cols = [c for c in enrich_df.columns if "__x__" in c]
mean_enrich = enrich_df[interaction_cols].mean()
mean_matrix = np.zeros((L, L))
for a, la in enumerate(LABELS):
    for b, lb in enumerate(LABELS):
        col = f"{la}__x__{lb}"
        if col in mean_enrich:
            mean_matrix[a, b] = mean_enrich[col]

mean_df = pd.DataFrame(mean_matrix, index=LABELS, columns=LABELS)
mean_df.to_csv(OUT_TAB / "mean_interaction_matrix.csv")

# %% ─── Figure 1: Mean interaction heatmap ────────────────────────────────────

fig, ax = plt.subplots(figsize=(13, 11))
sns.heatmap(mean_df, cmap="RdBu_r", center=0, vmin=-2, vmax=2, ax=ax,
            linewidths=0.3, cbar_kws={"label": "log2(observed/expected) contact enrichment"},
            annot=True, fmt=".2f", annot_kws={"size": 6})
ax.set_title("Mean pairwise cell-cell contact enrichment\n(radius-32, permutation-normalized, all tumor ROIs)")
ax.set_xlabel("Target cell type")
ax.set_ylabel("Source cell type")
plt.xticks(rotation=45, ha="right", fontsize=7)
plt.yticks(fontsize=7)
plt.tight_layout()
fig.savefig(OUT_FIG / "mean_interaction_heatmap.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Differential interactions: high vs low Gleason ───────────────────────

print("Testing differential interactions by clinical group...")

diff_results = []
for col in interaction_cols:
    sub = enrich_clin.dropna(subset=[col, "gleason_grp"])
    vals = sub[col].values
    gg = sub["gleason_grp"].astype(float).values
    rho, pval = spearmanr(gg, vals)
    diff_results.append({
        "interaction": col, "clin_var": "gleason_grp",
        "effect_size": rho, "pval": pval,
    })
    # stromogenic
    sub2 = enrich_clin.dropna(subset=[col, "stromogenic_smc_loss_reactive_stroma_present"])
    g0 = sub2[sub2.stromogenic_smc_loss_reactive_stroma_present == "no"][col].values
    g1 = sub2[sub2.stromogenic_smc_loss_reactive_stroma_present == "yes"][col].values
    if len(g0) >= 3 and len(g1) >= 3:
        stat, pval2 = mannwhitneyu(g0, g1, alternative="two-sided")
        rbc = 1 - 2 * stat / (len(g0) * len(g1))
        diff_results.append({"interaction": col, "clin_var": "stromogenic", "effect_size": rbc, "pval": pval2})
    # inflammation
    sub3 = enrich_clin.dropna(subset=[col, "inflammation"])
    g0i = sub3[sub3.inflammation == "no"][col].values
    g1i = sub3[sub3.inflammation == "yes"][col].values
    if len(g0i) >= 3 and len(g1i) >= 3:
        stat3, pval3 = mannwhitneyu(g0i, g1i, alternative="two-sided")
        rbc3 = 1 - 2 * stat3 / (len(g0i) * len(g1i))
        diff_results.append({"interaction": col, "clin_var": "inflammation", "effect_size": rbc3, "pval": pval3})

diff_df = pd.DataFrame(diff_results)
diff_df["fdr"] = multipletests(diff_df["pval"].fillna(1), method="fdr_bh")[1]
diff_df = diff_df.sort_values("fdr")
diff_df.to_csv(OUT_TAB / "differential_interactions.csv", index=False)

print(f"\nTop 15 differential interactions (FDR < 0.05):")
top = diff_df[diff_df.fdr < 0.05].head(15)
print(top[["interaction","clin_var","effect_size","fdr"]].to_string(index=False))

# %% ─── Figure 2: Differential interaction matrix for Gleason ────────────────

# Reshape Gleason Spearman rho into matrix
gleason_diff = diff_df[diff_df.clin_var == "gleason_grp"].set_index("interaction")["effect_size"]
gleason_matrix = np.zeros((L, L))
for a, la in enumerate(LABELS):
    for b, lb in enumerate(LABELS):
        col = f"{la}__x__{lb}"
        if col in gleason_diff:
            gleason_matrix[a, b] = gleason_diff[col]

gleason_df = pd.DataFrame(gleason_matrix, index=LABELS, columns=LABELS)

fig, ax = plt.subplots(figsize=(13, 11))
sns.heatmap(gleason_df, cmap="RdBu_r", center=0, vmin=-0.3, vmax=0.3, ax=ax,
            linewidths=0.3, cbar_kws={"label": "Spearman rho (correlation with Gleason grade)"},
            annot=True, fmt=".2f", annot_kws={"size": 6})
ax.set_title("Interaction enrichment correlation with Gleason grade group\n(Spearman rho; red = increases with grade)")
ax.set_xlabel("Target cell type")
ax.set_ylabel("Source cell type")
plt.xticks(rotation=45, ha="right", fontsize=7)
plt.yticks(fontsize=7)
plt.tight_layout()
fig.savefig(OUT_FIG / "gleason_differential_interactions.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 3: Top differential interaction pairs (volcano) ───────────────

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for ax, clin_var in zip(axes, ["gleason_grp", "stromogenic", "inflammation"]):
    sub = diff_df[diff_df.clin_var == clin_var].copy()
    x = sub["effect_size"].values
    y = -np.log10(sub["fdr"].values + 1e-15)
    colors = ["red" if f < 0.05 else "gray" for f in sub["fdr"]]
    ax.scatter(x, y, c=colors, s=8, alpha=0.6)
    ax.axhline(-np.log10(0.05), color="red", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)
    # Label top hits
    top_hits = sub[sub.fdr < 0.05].nlargest(5, "effect_size").index.tolist() + \
               sub[sub.fdr < 0.05].nsmallest(5, "effect_size").index.tolist()
    for idx in top_hits:
        row = sub.loc[idx]
        short_name = row["interaction"].replace("__x__", "↔").replace("stromal-","s-")\
                                        .replace("immune-","i-").replace("epithelial-","e-")\
                                        .replace("endothelial-","en-")
        ax.annotate(short_name, (row["effect_size"], -np.log10(row["fdr"] + 1e-15)),
                    fontsize=5, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("Effect size")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(f"Interaction differential — {clin_var}\n({(sub.fdr<0.05).sum()} significant)")

plt.tight_layout()
fig.savefig(OUT_FIG / "differential_interaction_volcanos.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Summary ───────────────────────────────────────────────────────────────

print(f"\n=== RQ6 Summary ===")
print(f"ROIs: {len(enrich_df)}, interaction pairs: {len(interaction_cols)}")
print(f"\nStrongest mean enrichments (self-contacts expected high):")
# Top off-diagonal enrichments
off_diag = mean_df.copy()
np.fill_diagonal(off_diag.values, np.nan)
stack = off_diag.stack().sort_values(ascending=False)
print(stack.head(10).to_string())
print(f"\nTop depleted off-diagonal:")
print(stack.tail(10).sort_values().to_string())
print(f"\nTotal significant differential interactions (FDR < 0.05): {(diff_df.fdr < 0.05).sum()}")
print(f"  vs gleason: {(diff_df[(diff_df.clin_var=='gleason_grp') & (diff_df.fdr<0.05)]).shape[0]}")
print(f"  vs stromogenic: {(diff_df[(diff_df.clin_var=='stromogenic') & (diff_df.fdr<0.05)]).shape[0]}")
print(f"  vs inflammation: {(diff_df[(diff_df.clin_var=='inflammation') & (diff_df.fdr<0.05)]).shape[0]}")
print(f"\nOutputs → {OUT_TAB}/ and {OUT_FIG}/")
