"""
RQ9: Spatial Niche Reconstruction via Local Composition K-Means

For each cell, computes a 17-dim local composition vector (meta-label fractions
in radius-32 neighborhood). K-means (k=18, matching the draft paper) identifies
recurrent spatial niches. Tests niche usage vs Gleason, stromogenic state, and
survival endpoints. Specifically tests whether the periglandular CAF1 niche
predicts clinical outcome, to reconcile with RQ4's counterintuitive direction.
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
from sklearn.cluster import MiniBatchKMeans
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ9")
OUT_TAB = Path("output/tables/RQ9")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32
K_NICHES = 18
MAX_CELLS_SAMPLE = 300_000   # cells for k-means fitting
CELLS_PER_ROI = 700          # max cells per ROI for fitting sample
RNG = np.random.default_rng(42)

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


def load_roi(sid: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (coords [N,2], labels [N]) for a single ROI."""
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    if "sample_id" in meta_df.index.names:
        meta_df = meta_df.xs(sid, level="sample_id") if sid in meta_df.index.get_level_values("sample_id") else meta_df
    spat_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")
    if "sample_id" in spat_df.index.names:
        spat_df = spat_df.xs(sid, level="sample_id") if sid in spat_df.index.get_level_values("sample_id") else spat_df
    # align
    idx = meta_df.index.intersection(spat_df.index)
    coords = spat_df.loc[idx, ["x", "y"]].to_numpy(dtype=np.float32)
    labels = meta_df.loc[idx, "meta_label"].to_numpy()
    return coords, labels


# Discover full label set from a subset of ROIs
all_labels: set = set()
for sid in sample_ids[:50]:
    _, lbl = load_roi(sid)
    all_labels.update(lbl)
LABELS = sorted(all_labels)
L = len(LABELS)
label_idx = {l: i for i, l in enumerate(LABELS)}
print(f"  {L} meta-labels: {LABELS}")


# %% ─── Local composition vectors ─────────────────────────────────────────────

def local_composition(coords: np.ndarray, labels: np.ndarray, radius: float) -> np.ndarray:
    """
    For each cell, compute the 17-dim meta-label fraction vector of its radius
    neighborhood (excluding self).

    Returns array of shape [N, L] with values summing to ≤1 (0 for isolated cells).
    """
    n_labels = L
    tree = cKDTree(coords)
    vecs = np.zeros((len(coords), n_labels), dtype=np.float32)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    if len(pairs) == 0:
        return vecs
    # accumulate both directions
    for a, b in pairs:
        vecs[a, label_idx[labels[b]]] += 1
        vecs[b, label_idx[labels[a]]] += 1
    # normalize
    totals = vecs.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1  # isolated cells stay zero
    vecs /= totals
    return vecs


print("Building local composition sample for k-means fitting...")
sample_vecs_list: list[np.ndarray] = []
sample_meta: list[dict] = []  # sid per sampled cell (for later assignment)

for sid in sample_ids:
    coords, labels = load_roi(sid)
    if len(coords) < 5:
        continue
    vecs = local_composition(coords, labels, RADIUS)

    # stratified subsample
    n_take = min(CELLS_PER_ROI, len(coords))
    idx = RNG.choice(len(coords), n_take, replace=False)
    sample_vecs_list.append(vecs[idx])

sample_vecs = np.vstack(sample_vecs_list)
print(f"  {len(sample_vecs):,} cells in fitting sample")

if len(sample_vecs) > MAX_CELLS_SAMPLE:
    keep = RNG.choice(len(sample_vecs), MAX_CELLS_SAMPLE, replace=False)
    sample_vecs = sample_vecs[keep]
    print(f"  Downsampled to {len(sample_vecs):,}")

# %% ─── K-means (k=18) ────────────────────────────────────────────────────────

print(f"Fitting MiniBatchKMeans k={K_NICHES}...")
km = MiniBatchKMeans(n_clusters=K_NICHES, random_state=42, batch_size=10_000, n_init=5)
km.fit(sample_vecs)
print("  Done.")

# Niche centers: [k, L]  (mean local composition per niche)
centers = pd.DataFrame(km.cluster_centers_, columns=LABELS)
centers.index.name = "niche"
centers.to_csv(OUT_TAB / "niche_centers.csv")
print(f"  Saved niche centers → {OUT_TAB / 'niche_centers.csv'}")

# %% ─── Per-ROI niche usage ────────────────────────────────────────────────────

print("Assigning all cells to niches and computing per-ROI usage...")
roi_niche_usage: list[dict] = []

for sid in sample_ids:
    coords, labels = load_roi(sid)
    if len(coords) < 5:
        continue
    vecs = local_composition(coords, labels, RADIUS)
    niche_assign = km.predict(vecs)
    counts = np.bincount(niche_assign, minlength=K_NICHES)
    fracs = counts / len(coords)
    row = {"sample_id": sid}
    for k in range(K_NICHES):
        row[f"niche_{k}"] = fracs[k]
    roi_niche_usage.append(row)

roi_niche_df = pd.DataFrame(roi_niche_usage).set_index("sample_id")
roi_niche_df.to_csv(OUT_TAB / "roi_niche_usage.csv")
print(f"  {len(roi_niche_df)} ROIs. Saved → {OUT_TAB / 'roi_niche_usage.csv'}")

# %% ─── Characterize niches ────────────────────────────────────────────────────

# Dominant meta-label and top-3 per niche
niche_names = []
for k in range(K_NICHES):
    top = centers.loc[k].nlargest(3)
    dominant = top.index[0]
    top3 = ", ".join([f"{c[:18]}({v:.2f})" for c, v in top.items()])
    niche_names.append(f"N{k}:{dominant[:12]}")
    print(f"  Niche {k:2d}: {top3}")

# Heatmap of niche centers
fig, ax = plt.subplots(figsize=(14, 8))
sns.heatmap(
    centers.T,
    ax=ax,
    cmap="YlOrRd",
    xticklabels=[f"N{k}" for k in range(K_NICHES)],
    yticklabels=LABELS,
)
ax.set_title(f"Spatial niche composition (k={K_NICHES}, radius-32)")
ax.set_xlabel("Niche")
ax.set_ylabel("Meta-label")
plt.tight_layout()
fig.savefig(OUT_FIG / "niche_composition_heatmap.pdf", bbox_inches="tight")
plt.close()
print(f"  Saved niche heatmap → {OUT_FIG / 'niche_composition_heatmap.pdf'}")

# Mean niche usage across all ROIs
mean_usage = roi_niche_df.mean()
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(range(K_NICHES), mean_usage.values)
ax.set_xticks(range(K_NICHES))
ax.set_xticklabels([f"N{k}" for k in range(K_NICHES)], rotation=45)
ax.set_ylabel("Mean fraction of cells")
ax.set_title("Mean niche usage across all ROIs")
plt.tight_layout()
fig.savefig(OUT_FIG / "mean_niche_usage.pdf", bbox_inches="tight")
plt.close()

# %% ─── Clinical associations: ROI level ──────────────────────────────────────

print("Testing niche usage vs clinical variables...")
roi_clin = clinical.loc[clinical.index.isin(roi_niche_df.index)].copy()
niche_cols = [f"niche_{k}" for k in range(K_NICHES)]
niche_usage = roi_niche_df.loc[roi_clin.index, niche_cols]

results_roi: list[dict] = []

for col in niche_cols:
    vals = niche_usage[col].values

    # Gleason numeric (gleason_grp is the grade group 1-5)
    gl = roi_clin["gleason_grp"].replace("", np.nan).dropna()
    gl_num = pd.to_numeric(gl, errors="coerce").dropna()
    common = gl_num.index.intersection(niche_usage.index)
    if len(common) > 20:
        rho, p = spearmanr(niche_usage.loc[common, col], gl_num.loc[common])
    else:
        rho, p = np.nan, np.nan
    results_roi.append({"niche": col, "variable": "gleason_spearman_rho", "stat": rho, "pval": p})

    # stromogenic
    for var in ["stromogenic", "inflammation", "cribriform"]:
        if var not in roi_clin.columns:
            continue
        grps = {g: niche_usage.loc[roi_clin[var] == g, col].dropna() for g in roi_clin[var].dropna().unique()}
        grps = {g: v for g, v in grps.items() if len(v) >= 5}
        if len(grps) == 2:
            g1, g2 = list(grps.values())
            stat, p = mannwhitneyu(g1, g2, alternative="two-sided")
            n = len(g1) + len(g2)
            rbc = 1 - 2 * stat / (len(g1) * len(g2))
            results_roi.append({"niche": col, "variable": var, "stat": rbc, "pval": p, "n": n})

roi_res_df = pd.DataFrame(results_roi)
# FDR per variable
for var in roi_res_df["variable"].unique():
    mask = roi_res_df["variable"] == var
    pvals = roi_res_df.loc[mask, "pval"].fillna(1).values
    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    roi_res_df.loc[mask, "fdr"] = fdr

roi_res_df.to_csv(OUT_TAB / "niche_roi_clinical_associations.csv", index=False)
sig = roi_res_df[roi_res_df["fdr"] < 0.1].sort_values("fdr")
print(f"  {len(sig)} significant ROI-level associations (FDR<0.1)")
print(sig[["niche", "variable", "stat", "pval", "fdr"]].to_string(index=False))

# %% ─── Patient-level niche usage ─────────────────────────────────────────────

print("\nAggregating to patient level...")
roi_niche_df["pat_id"] = clinical.loc[roi_niche_df.index.intersection(clinical.index), "pat_id"]
pat_niche = roi_niche_df.groupby("pat_id")[niche_cols].mean()
print(f"  {len(pat_niche)} patients")

# patient clinical
pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[
        ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time",
         "gleason_grp"]
    ]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")
pat_data = pat_niche.join(pat_clin, how="inner")
print(f"  {len(pat_data)} patients with clinical data")

# %% ─── Cox PH per niche ───────────────────────────────────────────────────────

print("Running Cox PH per niche vs survival endpoints...")
cox_results: list[dict] = []

for endpoint, duration_col in [
    ("psa_recurrence", "psa_progr_time"),
    ("clinical_progression", "clinical_progr_time"),
]:
    event_col = "psa_progr" if endpoint == "psa_recurrence" else "clinical_progr"
    sub = pat_data[[duration_col, event_col] + niche_cols].dropna(subset=[duration_col, event_col])
    sub = sub[sub[duration_col] > 0]
    n_events = sub[event_col].sum()
    print(f"  {endpoint}: n={len(sub)}, events={int(n_events)}")

    for col in niche_cols:
        cph_df = sub[[duration_col, event_col, col]].dropna()
        cph_df = cph_df.rename(columns={duration_col: "T", event_col: "E", col: "score"})
        if cph_df["score"].std() < 1e-8:
            continue
        try:
            cph = CoxPHFitter()
            cph.fit(cph_df, duration_col="T", event_col="E", show_progress=False)
            s = cph.summary.loc["score"]
            cox_results.append({
                "niche": col,
                "endpoint": endpoint,
                "hr": np.exp(s["coef"]),
                "ci_low": np.exp(s["coef lower 95%"]),
                "ci_high": np.exp(s["coef upper 95%"]),
                "p": s["p"],
                "n": len(cph_df),
                "n_events": int(cph_df["E"].sum()),
            })
        except Exception:
            pass

cox_df = pd.DataFrame(cox_results)
if len(cox_df):
    for ep in cox_df["endpoint"].unique():
        mask = cox_df["endpoint"] == ep
        _, fdr, _, _ = multipletests(cox_df.loc[mask, "p"].fillna(1).values, method="fdr_bh")
        cox_df.loc[mask, "fdr"] = fdr
    cox_df.to_csv(OUT_TAB / "niche_cox_results.csv", index=False)
    sig_cox = cox_df[cox_df["fdr"] < 0.1].sort_values("fdr")
    print(f"  {len(sig_cox)} significant niche Cox results (FDR<0.1)")
    print(sig_cox[["niche", "endpoint", "hr", "ci_low", "ci_high", "p", "fdr"]].to_string(index=False))

# %% ─── Forest plot: Cox HR per niche ─────────────────────────────────────────

for ep in ["psa_recurrence", "clinical_progression"]:
    ep_df = cox_df[cox_df["endpoint"] == ep].copy()
    if ep_df.empty:
        continue
    ep_df = ep_df.sort_values("hr")
    fig, ax = plt.subplots(figsize=(7, 8))
    y = range(len(ep_df))
    ax.scatter(ep_df["hr"], y, color=["tomato" if hr > 1 else "steelblue" for hr in ep_df["hr"]], zorder=3)
    ax.hlines(y, ep_df["ci_low"], ep_df["ci_high"], color="gray", linewidth=1)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(ep_df["niche"].str.replace("niche_", "N"), fontsize=8)
    ax.set_xlabel("Hazard Ratio (95% CI)")
    ax.set_title(f"Cox PH per niche — {ep}")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"cox_forest_{ep}.pdf", bbox_inches="tight")
    plt.close()

# %% ─── Identify periglandular CAF1 niche & reconcile with RQ4 ────────────────

print("\n─── Niche identification: periglandular CAF1 ───")
caf1_plus_label = next((l for l in LABELS if "CAF1" in l and "CD105+" in l), None)
luminal_label = next((l for l in LABELS if "luminal" in l.lower()), None)
print(f"  CAF1+ label: {caf1_plus_label}")
print(f"  Luminal label: {luminal_label}")

if caf1_plus_label and luminal_label:
    # Niche with highest combined CAF1+ + luminal fraction
    peri_score = centers[caf1_plus_label] + centers[luminal_label]
    peri_niche = int(peri_score.idxmax())
    print(f"  Periglandular CAF1 niche candidate: N{peri_niche}")
    print(f"    CAF1+={centers.loc[peri_niche, caf1_plus_label]:.3f}, "
          f"luminal={centers.loc[peri_niche, luminal_label]:.3f}")
    print(f"    Top composition: {centers.loc[peri_niche].nlargest(5).to_dict()}")

    if len(cox_df):
        peri_col = f"niche_{peri_niche}"
        sub = cox_df[cox_df["niche"] == peri_col]
        if not sub.empty:
            print(f"\n  Cox results for periglandular niche ({peri_col}):")
            print(sub[["endpoint", "hr", "ci_low", "ci_high", "p", "fdr"]].to_string(index=False))
        else:
            print(f"  No Cox results found for {peri_col}")

# %% ─── KM curves for top significant niche ───────────────────────────────────

if len(cox_df) and len(sig_cox):
    top_row = sig_cox.iloc[0]
    niche_col = top_row["niche"]
    ep = top_row["endpoint"]
    duration_col = "psa_progr_time" if ep == "psa_recurrence" else "clin_progr_time"
    event_col = "psa_progr" if ep == "psa_recurrence" else "clin_progr"

    sub = pat_data[[duration_col, event_col, niche_col]].dropna()
    sub = sub[sub[duration_col] > 0]
    median_val = sub[niche_col].median()
    hi = sub[sub[niche_col] >= median_val]
    lo = sub[sub[niche_col] < median_val]

    lr = logrank_test(hi[duration_col], lo[duration_col], hi[event_col], lo[event_col])

    fig, ax = plt.subplots(figsize=(6, 4))
    KaplanMeierFitter().fit(hi[duration_col], hi[event_col], label=f"High {niche_col}").plot_survival_function(ax=ax)
    KaplanMeierFitter().fit(lo[duration_col], lo[event_col], label=f"Low {niche_col}").plot_survival_function(ax=ax)
    ax.set_title(f"KM — {ep}\n{niche_col} high vs low (p={lr.p_value:.3f})")
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"km_{ep}_{niche_col}.pdf", bbox_inches="tight")
    plt.close()

print("\nRQ9 complete. Outputs saved to", OUT_TAB, "and", OUT_FIG)
