"""
RQ10: B-Cell Spatial Clustering as an Improved TLS Metric

Builds radius-32 spatial graphs restricted to B cells per ROI, extracts connected
component sizes as a measure of B-cell aggregation / TLS-like structure, and tests
association with clinical variables and survival.

Replaces the uninformative RQ3 TLS score (fraction of B cells with ≥1 T cell neighbor).
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
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ10")
OUT_TAB = Path("output/tables/RQ10")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32
BCELL_LABEL = "immune-B-cells"
LARGE_CLUSTER_MIN = 5  # threshold to call a "large" cluster

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


def connected_components_radius(coords: np.ndarray, radius: float) -> np.ndarray:
    """
    Union-Find connected components for points within `radius` of each other.
    Returns array of cluster labels (length = len(coords)).
    """
    n = len(coords)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    tree = cKDTree(coords)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    for a, b in pairs:
        union(int(a), int(b))

    return np.array([find(i) for i in range(n)])


# %% ─── Per-ROI B-cell clustering metrics ─────────────────────────────────────

print("Computing B-cell cluster metrics per ROI...")
roi_metrics: list[dict] = []

for sid in sample_ids:
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    spat_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")

    if "sample_id" in meta_df.index.names:
        meta_df = meta_df.xs(sid, level="sample_id") if sid in meta_df.index.get_level_values("sample_id") else meta_df
    if "sample_id" in spat_df.index.names:
        spat_df = spat_df.xs(sid, level="sample_id") if sid in spat_df.index.get_level_values("sample_id") else spat_df

    idx = meta_df.index.intersection(spat_df.index)
    labels = meta_df.loc[idx, "meta_label"].to_numpy()
    coords_all = spat_df.loc[idx, ["x", "y"]].to_numpy(dtype=np.float32)

    b_mask = labels == BCELL_LABEL
    n_bcells = b_mask.sum()

    if n_bcells < 2:
        roi_metrics.append({
            "sample_id": sid,
            "n_bcells": int(n_bcells),
            "n_clusters": 0,
            "max_cluster_size": int(n_bcells),
            "mean_cluster_size": float(n_bcells),
            "frac_in_large_cluster": 0.0,
            "bb_contact_density": 0.0,
        })
        continue

    b_coords = coords_all[b_mask]

    # Connected components
    comp_labels = connected_components_radius(b_coords, RADIUS)
    unique, counts = np.unique(comp_labels, return_counts=True)
    n_clusters = len(unique)
    max_size = int(counts.max())
    mean_size = float(counts.mean())
    frac_large = float((counts[counts >= LARGE_CLUSTER_MIN]).sum() / n_bcells)

    # B-B contact density: edges per B cell
    tree = cKDTree(b_coords)
    pairs = tree.query_pairs(RADIUS, output_type="ndarray")
    bb_density = len(pairs) / n_bcells if n_bcells > 0 else 0.0

    roi_metrics.append({
        "sample_id": sid,
        "n_bcells": int(n_bcells),
        "n_clusters": int(n_clusters),
        "max_cluster_size": max_size,
        "mean_cluster_size": mean_size,
        "frac_in_large_cluster": frac_large,
        "bb_contact_density": float(bb_density),
    })

roi_df = pd.DataFrame(roi_metrics).set_index("sample_id")
roi_df.to_csv(OUT_TAB / "roi_bcell_cluster_metrics.csv")
print(f"  {len(roi_df)} ROIs, {(roi_df['n_bcells'] > 0).sum()} with ≥1 B cell")
print(f"  {(roi_df['n_bcells'] >= 2).sum()} ROIs with ≥2 B cells (included in analysis)")
print(f"  Mean max cluster size: {roi_df['max_cluster_size'].mean():.1f} ± {roi_df['max_cluster_size'].std():.1f}")
print(f"  Mean B-B contact density: {roi_df['bb_contact_density'].mean():.2f}")

# %% ─── Clinical associations ─────────────────────────────────────────────────

print("\nTesting B-cell cluster metrics vs clinical variables...")
roi_clin = clinical.loc[clinical.index.isin(roi_df.index)].copy()
metric_cols = ["max_cluster_size", "mean_cluster_size", "frac_in_large_cluster", "bb_contact_density"]

results: list[dict] = []

for col in metric_cols:
    vals = roi_df.loc[roi_clin.index, col]

    # Gleason
    gl_num = pd.to_numeric(roi_clin["gleason_grp"].replace("", np.nan), errors="coerce").dropna()
    common = gl_num.index.intersection(vals.dropna().index)
    if len(common) > 20:
        rho, p = spearmanr(vals.loc[common], gl_num.loc[common])
        results.append({"metric": col, "variable": "gleason_spearman", "stat": rho, "pval": p, "n": len(common)})

    for var in ["inflammation", "stromogenic", "cribriform"]:
        if var not in roi_clin.columns:
            continue
        grps = {g: vals.loc[roi_clin[var] == g].dropna() for g in roi_clin[var].dropna().unique()}
        grps = {g: v for g, v in grps.items() if len(v) >= 5}
        if len(grps) == 2:
            g1, g2 = list(grps.values())
            stat, p = mannwhitneyu(g1, g2, alternative="two-sided")
            rbc = 1 - 2 * stat / (len(g1) * len(g2))
            results.append({"metric": col, "variable": var, "stat": rbc, "pval": p, "n": len(g1) + len(g2)})

res_df = pd.DataFrame(results)
_, fdr, _, _ = multipletests(res_df["pval"].fillna(1).values, method="fdr_bh")
res_df["fdr"] = fdr
res_df.to_csv(OUT_TAB / "bcell_cluster_roi_associations.csv", index=False)

sig = res_df[res_df["fdr"] < 0.1].sort_values("fdr")
print(f"  {len(sig)} significant associations (FDR<0.1):")
print(sig[["metric", "variable", "stat", "pval", "fdr", "n"]].to_string(index=False))

# %% ─── Boxplots: key associations ───────────────────────────────────────────

for var in ["inflammation", "stromogenic"]:
    if var not in roi_clin.columns:
        continue
    fig, axes = plt.subplots(1, len(metric_cols), figsize=(14, 4))
    for ax, col in zip(axes, metric_cols):
        plot_data = []
        labels_plot = []
        for g in sorted(roi_clin[var].dropna().unique()):
            v = roi_df.loc[roi_clin[roi_clin[var] == g].index, col].dropna()
            if len(v) >= 5:
                plot_data.append(v.values)
                labels_plot.append(str(g))
        if len(plot_data) == 2:
            ax.boxplot(plot_data, labels=labels_plot, showfliers=False)
            ax.set_title(col, fontsize=8)
            ax.set_ylabel(col, fontsize=7)
    fig.suptitle(f"B-cell cluster metrics vs {var}")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"bcell_clusters_vs_{var}.pdf", bbox_inches="tight")
    plt.close()

# %% ─── Patient-level aggregation ─────────────────────────────────────────────

print("\nAggregating to patient level...")
roi_df["pat_id"] = clinical.loc[roi_df.index.intersection(clinical.index), "pat_id"]

# Use max over ROIs (captures worst-case TLS-like infiltration)
pat_metrics = roi_df.groupby("pat_id")[metric_cols].max()
print(f"  {len(pat_metrics)} patients")

pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time"]]
    .first()
)
pat_data = pat_metrics.join(pat_clin, how="inner")
print(f"  {len(pat_data)} patients with clinical data")

# %% ─── Cox PH per metric ─────────────────────────────────────────────────────

print("Running Cox PH per B-cell cluster metric...")
cox_results: list[dict] = []

for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    sub = pat_data[[duration_col, event_col] + metric_cols].dropna(subset=[duration_col, event_col])
    sub = sub[sub[duration_col] > 0]
    print(f"  {endpoint}: n={len(sub)}, events={int(sub[event_col].sum())}")

    for col in metric_cols:
        cph_df = sub[[duration_col, event_col, col]].dropna()
        cph_df = cph_df.rename(columns={duration_col: "T", event_col: "E", col: "score"})
        if cph_df["score"].std() < 1e-8:
            continue
        try:
            cph = CoxPHFitter()
            cph.fit(cph_df, duration_col="T", event_col="E", show_progress=False)
            s = cph.summary.loc["score"]
            cox_results.append({
                "metric": col,
                "endpoint": endpoint,
                "hr": np.exp(s["coef"]),
                "ci_low": np.exp(s["coef lower 95%"]),
                "ci_high": np.exp(s["coef upper 95%"]),
                "p": s["p"],
                "n": len(cph_df),
                "n_events": int(cph_df["E"].sum()),
            })
        except Exception as e:
            print(f"    Cox failed for {col}/{endpoint}: {e}")

cox_df = pd.DataFrame(cox_results)
if len(cox_df):
    for ep in cox_df["endpoint"].unique():
        mask = cox_df["endpoint"] == ep
        _, fdr, _, _ = multipletests(cox_df.loc[mask, "p"].fillna(1).values, method="fdr_bh")
        cox_df.loc[mask, "fdr"] = fdr
    cox_df.to_csv(OUT_TAB / "bcell_cluster_cox.csv", index=False)
    print("\n  Cox PH results:")
    print(cox_df[["metric", "endpoint", "hr", "ci_low", "ci_high", "p", "fdr"]].sort_values("p").to_string(index=False))
    sig_cox = cox_df[cox_df["p"] < 0.1]
    if len(sig_cox):
        print(f"\n  Nominally significant (p<0.1): {len(sig_cox)}")

# %% ─── Compare with RQ3 TLS score ────────────────────────────────────────────

rq3_path = Path("output/tables/RQ3/roi_tls_scores.csv")
if rq3_path.exists():
    print("\nComparing with RQ3 TLS score...")
    tls_scores = pd.read_csv(rq3_path, index_col=0)
    tls_scores.index = tls_scores.index.astype(str)

    # correlation between TLS score and B-B contact density
    common = tls_scores.index.intersection(roi_df.index)
    if "tls_score" in tls_scores.columns and len(common) > 20:
        rho, p = spearmanr(tls_scores.loc[common, "tls_score"], roi_df.loc[common, "bb_contact_density"])
        print(f"  TLS score vs B-B contact density: Spearman rho={rho:.3f}, p={p:.3e}")
    if "bb_neighbor_frac" in tls_scores.columns and len(common) > 20:
        rho, p = spearmanr(tls_scores.loc[common, "bb_neighbor_frac"], roi_df.loc[common, "bb_contact_density"])
        print(f"  B-B neighbor fraction vs B-B contact density: Spearman rho={rho:.3f}, p={p:.3e}")

# %% ─── KM plot for best metric ───────────────────────────────────────────────

if len(cox_df):
    best = cox_df.sort_values("p").iloc[0]
    col, ep = best["metric"], best["endpoint"]
    duration_col = "psa_progr_time" if ep == "psa_recurrence" else "clinical_progr_time"
    event_col    = "psa_progr" if ep == "psa_recurrence" else "clinical_progr"

    sub = pat_data[[duration_col, event_col, col]].dropna()
    sub = sub[sub[duration_col] > 0]
    med = sub[col].median()
    hi = sub[sub[col] >= med]
    lo = sub[sub[col] < med]
    lr = logrank_test(hi[duration_col], lo[duration_col], hi[event_col], lo[event_col])

    fig, ax = plt.subplots(figsize=(6, 4))
    KaplanMeierFitter().fit(hi[duration_col], hi[event_col], label=f"High {col} (n={len(hi)})").plot_survival_function(ax=ax)
    KaplanMeierFitter().fit(lo[duration_col], lo[event_col], label=f"Low {col} (n={len(lo)})").plot_survival_function(ax=ax)
    ax.set_title(f"KM — {ep}\n{col} (log-rank p={lr.p_value:.3f}, HR={best['hr']:.2f})")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"km_{ep}_{col}.pdf", bbox_inches="tight")
    plt.close()
    print(f"\n  Saved KM plot → {OUT_FIG / f'km_{ep}_{col}.pdf'}")

print("\nRQ10 complete. Outputs saved to", OUT_TAB, "and", OUT_FIG)
