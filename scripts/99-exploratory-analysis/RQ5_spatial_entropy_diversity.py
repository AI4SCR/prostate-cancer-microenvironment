"""
RQ5: Spatial Entropy and Cell-Type Diversity

Computes global (ROI-level) and local (cell-level, radius-32 neighborhood) Shannon diversity
of cell types. Tests associations with clinical variables and survival.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr, kruskal
from statsmodels.stats.multitest import multipletests
from scipy.spatial import cKDTree
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ5")
OUT_TAB = Path("output/tables/RQ5")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32  # pixels

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

# %% ─── Helper: Shannon entropy ───────────────────────────────────────────────

def shannon(counts: np.ndarray) -> float:
    """Shannon entropy H = -sum(p * log(p)), base e."""
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0.0
    p = counts / counts.sum()
    return float(-np.sum(p * np.log(p)))


def local_shannon_per_cell(coords: np.ndarray, labels: np.ndarray,
                             radius: float, label_set: list) -> np.ndarray:
    """
    For each cell: compute Shannon entropy of cell type labels in its radius neighborhood
    (excluding self).
    Returns array of length n_cells.
    """
    tree = cKDTree(coords)
    n = len(coords)
    label_to_idx = {l: i for i, l in enumerate(label_set)}
    n_types = len(label_set)

    h_local = np.zeros(n)
    for i in range(n):
        neighbor_idxs = tree.query_ball_point(coords[i], r=radius)
        neighbor_idxs = [j for j in neighbor_idxs if j != i]
        if len(neighbor_idxs) == 0:
            h_local[i] = 0.0
            continue
        neighbor_labels = labels[neighbor_idxs]
        counts = np.zeros(n_types)
        for lbl in neighbor_labels:
            if lbl in label_to_idx:
                counts[label_to_idx[lbl]] += 1
        h_local[i] = shannon(counts)
    return h_local


# %% ─── Per-ROI global and local Shannon ──────────────────────────────────────

print("Computing Shannon diversity per ROI...")

# Load all cell type labels to get the full label set
all_labels_set: set = set()
for sid in sample_ids[:20]:  # sample to get label set
    m = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    m = m.loc[sid] if sid in m.index.get_level_values("sample_id") else m
    all_labels_set.update(m["meta_label"].unique())
label_set = sorted(all_labels_set)

roi_diversity = []
for i, sid in enumerate(sample_ids):
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{len(sample_ids)}")

    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    spatial_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")
    meta_df = meta_df.loc[sid] if sid in meta_df.index.get_level_values("sample_id") else meta_df
    spatial_df = spatial_df.loc[sid] if sid in spatial_df.index.get_level_values("sample_id") else spatial_df

    common = meta_df.index.intersection(spatial_df.index)
    if len(common) < 10:
        continue
    meta_df = meta_df.loc[common]
    spatial_df = spatial_df.loc[common]

    labels = meta_df["meta_label"].values
    coords = spatial_df[["x", "y"]].values
    n_cells = len(labels)

    # Global Shannon (ROI-level)
    label_arr = pd.Series(labels)
    counts = label_arr.value_counts().values
    h_global = shannon(counts)
    n_types = (counts > 0).sum()

    # Local Shannon (per-cell, radius-32)
    h_local = local_shannon_per_cell(coords, labels, RADIUS, label_set)
    h_local_mean = h_local.mean()
    h_local_median = float(np.median(h_local))
    h_local_sd = h_local.std()

    # Mean local Shannon per main cell type group
    main_labels = meta_df["main_group"].values
    h_by_group = {}
    for grp in ["epithelial", "stromal", "immune", "endothelial"]:
        mask = main_labels == grp
        if mask.sum() >= 3:
            h_by_group[f"h_local_mean_{grp}"] = h_local[mask].mean()

    record = {
        "sample_id": sid,
        "n_cells": n_cells,
        "n_types_observed": n_types,
        "h_global": h_global,
        "h_local_mean": h_local_mean,
        "h_local_median": h_local_median,
        "h_local_sd": h_local_sd,
    }
    record.update(h_by_group)
    roi_diversity.append(record)

diversity_df = pd.DataFrame(roi_diversity).set_index("sample_id")
print(f"  Done: {len(diversity_df)} ROIs")

# Merge with clinical
div_clin = diversity_df.join(clinical[["pat_id", "gleason_grp", "gleason_pattern_tma_core",
                                        "psa_progr", "psa_progr_time", "clinical_progr",
                                        "clinical_progr_time", "stromogenic_smc_loss_reactive_stroma_present",
                                        "inflammation", "cribriform"]])
div_clin.to_csv(OUT_TAB / "roi_shannon_diversity.csv")

# %% ─── Association tests ─────────────────────────────────────────────────────

print("\nRunning association tests...")

DIVERSITY_METRICS = [
    "h_global", "h_local_mean", "h_local_median", "h_local_sd",
    "h_local_mean_epithelial", "h_local_mean_stromal", "h_local_mean_immune",
]
CLINICAL_VARS = {
    "gleason_grp": "spearman",
    "stromogenic_smc_loss_reactive_stroma_present": "mannwhitney",
    "inflammation": "mannwhitney",
    "cribriform": "mannwhitney",
    "psa_progr": "mannwhitney",
    "clinical_progr": "mannwhitney",
}

results = []
for clin_col, test in CLINICAL_VARS.items():
    for metric in DIVERSITY_METRICS:
        if metric not in div_clin.columns:
            continue
        sub = div_clin.dropna(subset=[metric, clin_col])
        vals = sub[metric].values
        grp = sub[clin_col].values
        if test == "spearman":
            rho, pval = spearmanr(grp.astype(float), vals)
            es = rho
        else:
            g0 = vals[grp == (0 if grp.dtype in [float, int] else "no")]
            g1 = vals[grp == (1 if grp.dtype in [float, int] else "yes")]
            if len(g0) < 3 or len(g1) < 3:
                continue
            stat, pval = mannwhitneyu(g0, g1, alternative="two-sided")
            es = 1 - 2 * stat / (len(g0) * len(g1))
        results.append({"metric": metric, "clin_var": clin_col, "effect_size": es, "pval": pval})

res_df = pd.DataFrame(results)
res_df["fdr"] = multipletests(res_df["pval"].fillna(1), method="fdr_bh")[1]
res_df = res_df.sort_values("fdr")
res_df.to_csv(OUT_TAB / "shannon_association_results.csv", index=False)

print("Significant associations (FDR < 0.1):")
print(res_df[res_df.fdr < 0.1].to_string(index=False))

# %% ─── Figure 1: Global vs local diversity by Gleason grade ─────────────────

gleason_map = {1.0: "GG1", 2.0: "GG2", 3.0: "GG3", 4.0: "GG4", 5.0: "GG5"}
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, (metric, label) in zip(axes, [
    ("h_global", "Global Shannon diversity (ROI)"),
    ("h_local_mean", "Mean local Shannon diversity (cell, radius-32)"),
    ("h_local_sd", "SD local Shannon diversity (spatial heterogeneity)"),
]):
    sub = div_clin.dropna(subset=[metric, "gleason_grp"]).copy()
    sub["gg"] = sub["gleason_grp"].map(gleason_map)
    rho, pval = spearmanr(sub["gleason_grp"].astype(float).values, sub[metric].values)
    sns.boxplot(data=sub, x="gg", y=metric, ax=ax, palette="Blues",
                order=["GG1", "GG2", "GG3", "GG4", "GG5"])
    ax.set_title(f"{label}\nSpearman rho={rho:.3f}, p={pval:.3e}")
    ax.set_xlabel("Gleason grade group")
    ax.set_ylabel(label)

plt.tight_layout()
fig.savefig(OUT_FIG / "shannon_by_gleason.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 2: Diversity by stromogenic, inflammation, cribriform ──────────

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()

plot_combos = [
    ("h_local_mean", "stromogenic_smc_loss_reactive_stroma_present", "Stromogenic", ["no","yes"]),
    ("h_local_mean", "inflammation", "Inflammation", ["no","yes"]),
    ("h_local_mean", "cribriform", "Cribriform", ["no","yes"]),
    ("h_global",     "stromogenic_smc_loss_reactive_stroma_present", "Stromogenic", ["no","yes"]),
    ("h_global",     "inflammation", "Inflammation", ["no","yes"]),
    ("h_global",     "cribriform", "Cribriform", ["no","yes"]),
]

for ax, (metric, col, col_label, order) in zip(axes, plot_combos):
    sub = div_clin.dropna(subset=[metric, col])
    g0 = sub[sub[col] == order[0]][metric].values
    g1 = sub[sub[col] == order[1]][metric].values
    if len(g0) > 3 and len(g1) > 3:
        stat, pval = mannwhitneyu(g0, g1, alternative="two-sided")
        rbc = 1 - 2*stat/(len(g0)*len(g1))
    else:
        pval, rbc = np.nan, np.nan
    sub2 = sub.copy()
    sns.boxplot(data=sub2, x=col, y=metric, order=order, ax=ax,
                palette={"no": "#aec6cf", "yes": "#ff7f7f"})
    ax.set_title(f"{'Local' if 'local' in metric else 'Global'} Shannon vs {col_label}\np={pval:.3f}, rbc={rbc:.3f}" if not np.isnan(pval) else "")
    ax.set_xlabel(col_label)
    ax.set_ylabel(metric)

plt.tight_layout()
fig.savefig(OUT_FIG / "shannon_by_pathology.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 3: Local diversity per cell type group ────────────────────────

group_metrics = [c for c in diversity_df.columns if c.startswith("h_local_mean_")]
if group_metrics:
    melt_data = []
    for sid, row in div_clin.iterrows():
        for gm in group_metrics:
            if pd.notna(row.get(gm)):
                grp = gm.replace("h_local_mean_", "")
                gg = gleason_map.get(row.get("gleason_grp"), "unknown")
                melt_data.append({"group": grp, "h_local": row[gm], "gleason": gg})
    melt_df = pd.DataFrame(melt_data)

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(data=melt_df, x="group", y="h_local", hue="gleason",
                palette="Blues", hue_order=["GG1","GG2","GG3","GG4","GG5"], ax=ax)
    ax.set_title("Local Shannon diversity by cell type group and Gleason grade")
    ax.set_xlabel("Cell type group")
    ax.set_ylabel("Mean local Shannon diversity")
    ax.legend(title="Gleason", bbox_to_anchor=(1.05, 1), fontsize=8)
    plt.tight_layout()
    fig.savefig(OUT_FIG / "shannon_by_group_and_gleason.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

# %% ─── Figure 4: Global vs local Shannon (are they correlated?) ─────────────

fig, ax = plt.subplots(figsize=(7, 5))
rho, pval = spearmanr(div_clin["h_global"].dropna(), div_clin["h_local_mean"].dropna(),
                       nan_policy="omit")
ax.scatter(div_clin["h_global"], div_clin["h_local_mean"], alpha=0.3, s=15,
           c=div_clin["gleason_grp"], cmap="RdYlBu_r")
ax.set_xlabel("Global Shannon diversity (ROI)")
ax.set_ylabel("Mean local Shannon diversity (cell)")
ax.set_title(f"Global vs local diversity\nSpearman rho={rho:.3f}, p={pval:.3e}")
plt.tight_layout()
fig.savefig(OUT_FIG / "global_vs_local_shannon.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Patient-level survival ────────────────────────────────────────────────

print("\nPatient-level survival analysis...")

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

pat_clinical = clinical.groupby("pat_id").first()[
    ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time"]
]

# Max-pool diversity to patient level
div_pat = div_clin[["pat_id"] + DIVERSITY_METRICS].dropna(subset=["pat_id"])
div_pat_agg = div_pat.groupby("pat_id")[DIVERSITY_METRICS].mean()
div_pat_all = div_pat_agg.join(pat_clinical)

cox_results = []
for metric in ["h_global", "h_local_mean", "h_local_sd"]:
    if metric not in div_pat_all.columns:
        continue
    for event_col, time_col, label in [
        ("psa_progr", "psa_progr_time", "PSA recurrence"),
        ("clinical_progr", "clinical_progr_time", "Clinical progression"),
    ]:
        sub = div_pat_all.dropna(subset=[metric, event_col, time_col])
        if len(sub) < 20 or sub[event_col].sum() < 5:
            continue
        try:
            cph = CoxPHFitter()
            df = sub[[metric, event_col, time_col]].copy()
            df[metric] = (df[metric] - df[metric].mean()) / (df[metric].std() + 1e-10)
            cph.fit(df, duration_col=time_col, event_col=event_col, show_progress=False)
            r = cph.summary.loc[metric]
            cox_results.append({
                "metric": metric, "endpoint": label,
                "HR": float(np.exp(r["coef"])),
                "HR_lo": float(np.exp(r["coef lower 95%"])),
                "HR_hi": float(np.exp(r["coef upper 95%"])),
                "p": float(r["p"]), "n": len(sub), "n_events": int(sub[event_col].sum()),
            })
        except Exception as e:
            print(f"  Cox failed {metric}/{label}: {e}")

if cox_results:
    cox_df = pd.DataFrame(cox_results).sort_values("p")
    cox_df.to_csv(OUT_TAB / "shannon_cox_results.csv", index=False)
    print("\nCox PH results (Shannon diversity vs survival):")
    print(cox_df.to_string(index=False))

# %% ─── Summary ───────────────────────────────────────────────────────────────

print(f"\n=== RQ5 Summary ===")
print(f"ROIs processed: {len(diversity_df)}")
print(f"\nMean global Shannon: {div_clin.h_global.mean():.3f}")
print(f"Mean local Shannon:  {div_clin.h_local_mean.mean():.3f}")
print(f"\nSignificant associations (FDR < 0.1): {(res_df.fdr < 0.1).sum()}")
print(res_df[res_df.fdr < 0.1][["metric","clin_var","effect_size","fdr"]].to_string(index=False))
print(f"\nOutputs → {OUT_TAB}/ and {OUT_FIG}/")
