"""
RQ18: Dense Periglandular CAF1 Niche vs Diffuse CAF1 Proximity

Resolves the counterintuitive Finding 3 (more periglandular CAF1 → better outcome).
Distinguishes two modes:
  (a) Dense/organized: CAF1-CD105+ cells forming tight clusters around gland perimeters
  (b) Diffuse/scattered: CAF1-CD105+ cells present near glands but loosely distributed

Uses DBSCAN clustering of CAF1-CD105+ cells within 32px of luminal epithelium to
identify dense periglandular clusters, then computes:
  - Fraction of periglandular CAF1 in dense clusters (≥5 cells)
  - Max periglandular cluster size per ROI
  - "Niche score": proportion of CAF1-CD105+ cells in dense periglandular clusters

Tests these niche-level metrics vs the continuous proximity fraction (RQ4),
and vs clinical outcome.
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
from sklearn.cluster import DBSCAN
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ18")
OUT_TAB = Path("output/tables/RQ18")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32
# Dense cluster: DBSCAN epsilon and min_samples
DBSCAN_EPS = 40       # pixels — cells within 40px considered neighbors in cluster
DBSCAN_MIN = 4        # minimum 4 CAF1 cells to form a dense cluster
DENSE_MIN_SIZE = 5    # minimum cluster size to count as a "dense niche"

CAF1_PLUS_LABEL = "stromal-CAF1(CD105+)"   # fine-grained label
LUMINAL_LABELS  = {"epithelial-luminal", "epithelial-luminal(Ki67+)"}
TUMOR_LABELS    = {"epithelial-(ERG+CD44+)", "epithelial-luminal(ERG+)",
                   "epithelial-luminal(ERG+p53+)"}
ALL_GLAND_LABELS = LUMINAL_LABELS | TUMOR_LABELS

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
print(f"  {len(sample_ids)} tumor ROIs")

# %% ─── Per-ROI dense periglandular CAF1 niche metrics ────────────────────────

print("Computing dense periglandular CAF1 niche metrics...")
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

    caf1_mask  = labels == CAF1_PLUS_LABEL
    gland_mask = np.isin(labels, list(ALL_GLAND_LABELS))
    n_caf1     = caf1_mask.sum()
    n_gland    = gland_mask.sum()

    row: dict = {"sample_id": sid, "n_caf1": int(n_caf1), "n_gland": int(n_gland)}

    if n_caf1 < 2 or n_gland < 2:
        row.update({
            "n_periglandular_caf1": 0,
            "periglandular_frac": 0.0,
            "n_dense_clusters": 0,
            "max_dense_cluster_size": 0,
            "frac_caf1_in_dense_niche": 0.0,
            "dense_niche_score": 0.0,
        })
        roi_records.append(row)
        continue

    caf1_coords = coords[caf1_mask]
    gland_coords = coords[gland_mask]

    # Periglandular: CAF1+ cells within RADIUS of any gland cell
    tree_gland = cKDTree(gland_coords)
    dists, _ = tree_gland.query(caf1_coords, k=1)
    peri_mask = dists <= RADIUS
    n_peri = int(peri_mask.sum())
    peri_frac = n_peri / n_caf1

    # Dense clustering of periglandular CAF1 cells via DBSCAN
    dense_cluster_sizes: list[int] = []
    n_in_dense = 0
    if n_peri >= DBSCAN_MIN:
        peri_coords = caf1_coords[peri_mask]
        db = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN).fit(peri_coords)
        cluster_ids = db.labels_
        for cl in set(cluster_ids):
            if cl == -1:
                continue
            size = int((cluster_ids == cl).sum())
            if size >= DENSE_MIN_SIZE:
                dense_cluster_sizes.append(size)
                n_in_dense += size

    n_dense = len(dense_cluster_sizes)
    max_dense = max(dense_cluster_sizes) if dense_cluster_sizes else 0
    frac_in_dense = n_in_dense / n_caf1  # fraction of ALL CAF1+ cells in a dense periglandular cluster

    # Dense niche score: combines presence and size
    dense_niche_score = max_dense / (n_caf1 + 1)

    row.update({
        "n_periglandular_caf1": n_peri,
        "periglandular_frac": float(peri_frac),
        "n_dense_clusters": n_dense,
        "max_dense_cluster_size": max_dense,
        "frac_caf1_in_dense_niche": float(frac_in_dense),
        "dense_niche_score": float(dense_niche_score),
    })
    roi_records.append(row)

roi_df = pd.DataFrame(roi_records).set_index("sample_id")
roi_df.to_csv(OUT_TAB / "roi_dense_caf1_niche.csv")
print(f"  {len(roi_df)} ROIs")
print(f"  ROIs with ≥1 dense CAF1 cluster: {(roi_df['n_dense_clusters'] > 0).sum()}")
print(f"  Mean periglandular fraction: {roi_df['periglandular_frac'].mean():.3f}")
print(f"  Mean frac in dense niche: {roi_df['frac_caf1_in_dense_niche'].mean():.3f}")

# Correlation: dense niche score vs continuous periglandular fraction (RQ4 metric)
rho, p = spearmanr(roi_df["periglandular_frac"], roi_df["frac_caf1_in_dense_niche"])
print(f"\n  periglandular_frac vs frac_in_dense_niche: rho={rho:.3f}, p={p:.3e}")
rho2, p2 = spearmanr(roi_df["periglandular_frac"], roi_df["max_dense_cluster_size"])
print(f"  periglandular_frac vs max_dense_cluster_size: rho={rho2:.3f}, p={p2:.3e}")

# %% ─── Clinical associations ─────────────────────────────────────────────────

print("\nTesting dense niche metrics vs clinical variables...")
roi_clin = clinical.loc[clinical.index.isin(roi_df.index)].copy()

metric_cols = ["periglandular_frac", "frac_caf1_in_dense_niche",
               "max_dense_cluster_size", "dense_niche_score"]
results: list[dict] = []

for col in metric_cols:
    vals = roi_df.loc[roi_clin.index, col].dropna()
    gl = pd.to_numeric(roi_clin.loc[vals.index, "gleason_grp"].replace("", np.nan),
                       errors="coerce").dropna()
    common = vals.index.intersection(gl.index)
    if len(common) > 20:
        rho, p = spearmanr(vals.loc[common], gl.loc[common])
        results.append({"metric": col, "variable": "gleason", "stat": rho, "pval": p})

    for var in ["inflammation", "stromogenic", "cribriform"]:
        if var not in roi_clin.columns:
            continue
        grps = {g: vals.loc[roi_clin.loc[vals.index, var] == g].dropna()
                for g in roi_clin[var].dropna().unique()}
        grps = {g: v for g, v in grps.items() if len(v) >= 5}
        if len(grps) == 2:
            g1, g2 = list(grps.values())
            stat, p = mannwhitneyu(g1, g2, alternative="two-sided")
            rbc = 1 - 2 * stat / (len(g1) * len(g2))
            results.append({"metric": col, "variable": var, "stat": rbc, "pval": p})

res_df = pd.DataFrame(results)
_, fdr, _, _ = multipletests(res_df["pval"].fillna(1), method="fdr_bh")
res_df["fdr"] = fdr
res_df.to_csv(OUT_TAB / "dense_caf1_associations.csv", index=False)
print(res_df.sort_values("pval")[["metric","variable","stat","pval","fdr"]].to_string(index=False))

# %% ─── Patient-level Cox PH ──────────────────────────────────────────────────

print("\nPatient-level Cox PH...")
roi_df["pat_id"] = clinical.loc[roi_df.index.intersection(clinical.index), "pat_id"]
pat_metrics = roi_df.groupby("pat_id")[metric_cols].mean()
pat_metrics["frac_caf1_in_dense_niche_max"] = roi_df.groupby("pat_id")["frac_caf1_in_dense_niche"].max()
pat_metrics["dense_niche_score_max"] = roi_df.groupby("pat_id")["dense_niche_score"].max()

pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[["psa_progr", "psa_progr_time", "clinical_progr",
                         "clinical_progr_time", "gleason_grp"]]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")
pat_data = pat_metrics.join(pat_clin, how="inner")
all_cols = list(pat_metrics.columns)

cox_results: list[dict] = []
for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    sub = pat_data[[duration_col, event_col, "gleason_grp"] + all_cols].dropna(
        subset=[duration_col, event_col, "gleason_grp"]
    )
    sub = sub[sub[duration_col] > 0]
    print(f"  {endpoint}: n={len(sub)}, events={int(sub[event_col].sum())}")

    for col in all_cols:
        cph_df = sub[[duration_col, event_col, col]].dropna().rename(
            columns={duration_col: "T", event_col: "E", col: "score"}
        )
        if cph_df["score"].std() < 1e-8:
            continue
        try:
            cph_uni = CoxPHFitter()
            cph_uni.fit(cph_df, duration_col="T", event_col="E", show_progress=False)
            s_uni = cph_uni.summary.loc["score"]

            cph_adj_df = sub[[duration_col, event_col, "gleason_grp", col]].dropna().rename(
                columns={duration_col: "T", event_col: "E", col: "score"}
            )
            cph_adj = CoxPHFitter()
            cph_adj.fit(cph_adj_df, duration_col="T", event_col="E", show_progress=False)
            s_adj = cph_adj.summary.loc["score"]

            cox_results.append({
                "metric": col, "endpoint": endpoint,
                "hr_uni": np.exp(s_uni["coef"]), "p_uni": s_uni["p"],
                "hr_adj": np.exp(s_adj["coef"]),
                "ci_low_adj": np.exp(s_adj["coef lower 95%"]),
                "ci_high_adj": np.exp(s_adj["coef upper 95%"]),
                "p_adj": s_adj["p"],
                "n": len(cph_df),
            })
        except Exception:
            pass

cox_df = pd.DataFrame(cox_results)
if len(cox_df):
    cox_df.to_csv(OUT_TAB / "dense_caf1_cox.csv", index=False)
    print("\n  Cox results (sorted by Gleason-adjusted p):")
    print(cox_df.sort_values("p_adj")[
        ["metric", "endpoint", "hr_uni", "p_uni", "hr_adj", "ci_low_adj", "ci_high_adj", "p_adj"]
    ].to_string(index=False))

# %% ─── Figures ───────────────────────────────────────────────────────────────

# Scatter: periglandular_frac vs frac_in_dense_niche
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, color_col, title in zip(axes, ["gleason_grp", "cribriform"], ["Gleason GG", "Cribriform"]):
    sub = roi_clin.copy()
    sub["periglandular_frac"] = roi_df["periglandular_frac"]
    sub["frac_caf1_in_dense_niche"] = roi_df["frac_caf1_in_dense_niche"]
    sub = sub.dropna(subset=["periglandular_frac", "frac_caf1_in_dense_niche", color_col])
    if pd.api.types.is_numeric_dtype(sub[color_col]) or sub[color_col].apply(lambda x: str(x).replace('.','').isdigit()).all():
        sc = ax.scatter(sub["periglandular_frac"], sub["frac_caf1_in_dense_niche"],
                        c=pd.to_numeric(sub[color_col], errors="coerce"), cmap="RdYlBu_r",
                        s=10, alpha=0.6)
        plt.colorbar(sc, ax=ax, label=title)
    else:
        for cat in sorted(sub[color_col].unique()):
            m = sub[sub[color_col] == cat]
            ax.scatter(m["periglandular_frac"], m["frac_caf1_in_dense_niche"],
                       label=str(cat), s=10, alpha=0.6)
        ax.legend(fontsize=7, title=title)
    ax.set_xlabel("Periglandular fraction (RQ4 metric)")
    ax.set_ylabel("Fraction CAF1+ in dense niche (RQ18)")
    ax.set_title(f"Dense niche vs diffuse proximity — {title}")
plt.tight_layout()
fig.savefig(OUT_FIG / "dense_vs_diffuse_scatter.pdf", bbox_inches="tight")
plt.close()

# KM for best metric
if len(cox_df):
    best = cox_df.sort_values("p_adj").iloc[0]
    col, ep = best["metric"], best["endpoint"]
    duration_col = "psa_progr_time" if ep == "psa_recurrence" else "clinical_progr_time"
    event_col    = "psa_progr"       if ep == "psa_recurrence" else "clinical_progr"
    sub = pat_data[[duration_col, event_col, col]].dropna()
    sub = sub[sub[duration_col] > 0]
    med = sub[col].median()
    hi, lo = sub[sub[col] >= med], sub[sub[col] < med]
    lr = logrank_test(hi[duration_col], lo[duration_col], hi[event_col], lo[event_col])
    fig, ax = plt.subplots(figsize=(6, 4))
    KaplanMeierFitter().fit(hi[duration_col], hi[event_col],
                            label=f"High {col[:25]} (n={len(hi)})").plot_survival_function(ax=ax)
    KaplanMeierFitter().fit(lo[duration_col], lo[event_col],
                            label=f"Low (n={len(lo)})").plot_survival_function(ax=ax)
    ax.set_title(f"KM — {ep}\n{col} (log-rank p={lr.p_value:.3f})")
    ax.set_xlabel("Time (months)"); ax.set_ylabel("Survival probability")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"km_{ep}_{col[:30]}.pdf", bbox_inches="tight")
    plt.close()

print("\nRQ18 complete. Outputs →", OUT_TAB, "and", OUT_FIG)
