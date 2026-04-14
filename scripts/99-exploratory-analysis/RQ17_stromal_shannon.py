"""
RQ17: Stromal-Only Shannon Diversity — Disentangling Architecture from Inflammation

Tests whether spatial Shannon diversity (Finding 6: HR=1.38 for PSA recurrence) reflects
genuine tissue architectural disorganization or is simply driven by immune infiltration.

Approach:
  - Compute Shannon diversity restricted to non-immune cells (epithelial + stromal + endothelial)
  - Compare with full-composition Shannon (RQ5)
  - Test stromal-only Shannon vs survival (Gleason-adjusted)
  - If stromal Shannon is still prognostic → genuine architectural signal
  - If null → the RQ5 finding was driven by inflammation co-occurrence

Also computes epithelial-only and stromal-only spatial diversity separately.
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
OUT_FIG = Path("output/figures/RQ17")
OUT_TAB = Path("output/tables/RQ17")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32

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


def shannon(counts: np.ndarray) -> float:
    p = counts[counts > 0]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


def local_shannon_subset(
    coords: np.ndarray,
    meta_labels: np.ndarray,
    subset_mask: np.ndarray,
    radius: float,
    label_set: list[str],
) -> np.ndarray:
    """
    For each cell in subset_mask, compute Shannon entropy of meta_label distribution
    in its radius neighborhood (all cell types, not just subset).
    Returns array of length subset_mask.sum().
    """
    label_to_idx = {l: i for i, l in enumerate(label_set)}
    n_types = len(label_set)
    sub_coords = coords[subset_mask]
    tree = cKDTree(coords)
    h = np.zeros(len(sub_coords))
    for i, c in enumerate(sub_coords):
        nbrs = tree.query_ball_point(c, r=radius)
        nbrs = [j for j in nbrs if not subset_mask[j] or coords[j].tolist() != c.tolist()]
        if not nbrs:
            h[i] = 0.0
            continue
        cnts = np.zeros(n_types)
        for j in nbrs:
            lbl = meta_labels[j]
            if lbl in label_to_idx:
                cnts[label_to_idx[lbl]] += 1
        h[i] = shannon(cnts)
    return h


# %% ─── Per-ROI stromal-only and epithelial-only Shannon ─────────────────────

print("Computing cell-type-restricted Shannon diversity per ROI...")
roi_records: list[dict] = []

for sid in sample_ids:
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    spat_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")

    if "sample_id" in meta_df.index.names:
        meta_df = meta_df.xs(sid, level="sample_id") if sid in meta_df.index.get_level_values("sample_id") else meta_df
    if "sample_id" in spat_df.index.names:
        spat_df = spat_df.xs(sid, level="sample_id") if sid in spat_df.index.get_level_values("sample_id") else spat_df

    idx = meta_df.index.intersection(spat_df.index)
    meta_labels = meta_df.loc[idx, "meta_label"].to_numpy()
    main_groups = meta_df.loc[idx, "main_group"].to_numpy()
    coords = spat_df.loc[idx, ["x", "y"]].to_numpy(dtype=np.float32)

    label_set = sorted(set(meta_labels))

    # Global Shannon (all cells) — for comparison with RQ5
    counts_all = np.array([(meta_labels == l).sum() for l in label_set], dtype=float)
    h_global = shannon(counts_all)

    # Global Shannon restricted to non-immune cells
    non_immune_mask = main_groups != "immune"
    non_immune_labels = meta_labels[non_immune_mask]
    non_immune_label_set = sorted(set(non_immune_labels))
    if len(non_immune_label_set) >= 2:
        cnt_ni = np.array([(non_immune_labels == l).sum() for l in non_immune_label_set], dtype=float)
        h_non_immune = shannon(cnt_ni)
    else:
        h_non_immune = np.nan

    # Global stromal-only Shannon
    strom_mask  = main_groups == "stromal"
    strom_labels = meta_labels[strom_mask]
    if len(set(strom_labels)) >= 2:
        cnt_s = np.array([(strom_labels == l).sum() for l in set(strom_labels)], dtype=float)
        h_stromal = shannon(cnt_s)
    else:
        h_stromal = np.nan

    # Local Shannon (per non-immune cell) — neighborhood includes all cell types
    h_local_non_immune = np.nan
    h_local_stromal    = np.nan
    strom_idx = np.where(strom_mask)[0]
    non_immune_idx = np.where(non_immune_mask)[0]

    if len(non_immune_idx) >= 5 and len(label_set) >= 2:
        tree = cKDTree(coords)
        label_to_idx = {l: i for i, l in enumerate(label_set)}
        n_types = len(label_set)

        h_vals_ni = np.zeros(len(non_immune_idx))
        for k, ci in enumerate(non_immune_idx):
            nbrs = [j for j in tree.query_ball_point(coords[ci], r=RADIUS) if j != ci]
            if not nbrs:
                continue
            cnts = np.zeros(n_types)
            for j in nbrs:
                lbl = meta_labels[j]
                if lbl in label_to_idx:
                    cnts[label_to_idx[lbl]] += 1
            h_vals_ni[k] = shannon(cnts)
        h_local_non_immune = float(h_vals_ni.mean())

        if len(strom_idx) >= 3:
            h_vals_s = h_vals_ni[np.isin(non_immune_idx, strom_idx)]
            h_local_stromal = float(h_vals_s.mean()) if len(h_vals_s) > 0 else np.nan

    roi_records.append({
        "sample_id": sid,
        "h_global": h_global,
        "h_non_immune": h_non_immune,
        "h_stromal_global": h_stromal,
        "h_local_non_immune": h_local_non_immune,
        "h_local_stromal": h_local_stromal,
        "immune_frac": (main_groups == "immune").mean(),
        "n_cells": len(idx),
    })

roi_df = pd.DataFrame(roi_records).set_index("sample_id")
roi_df.to_csv(OUT_TAB / "roi_restricted_shannon.csv")
print(f"  {len(roi_df)} ROIs")

# Correlation between global and non-immune Shannon
from scipy.stats import spearmanr as _sp
rho, p = _sp(roi_df["h_global"].dropna(), roi_df["h_non_immune"].dropna(),
             nan_policy="omit")
print(f"  Global vs non-immune Shannon: rho={rho:.3f}, p={p:.2e}")
rho2, p2 = _sp(roi_df["h_global"].dropna(), roi_df["immune_frac"].dropna(), nan_policy="omit")
print(f"  Global Shannon vs immune fraction: rho={rho2:.3f}, p={p2:.2e}")

# %% ─── Clinical associations ─────────────────────────────────────────────────

print("\nTesting restricted Shannon metrics vs clinical variables...")
roi_clin = clinical.loc[clinical.index.isin(roi_df.index)].copy()
metric_cols = ["h_global", "h_non_immune", "h_stromal_global",
               "h_local_non_immune", "h_local_stromal"]

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
res_df.to_csv(OUT_TAB / "roi_restricted_shannon_associations.csv", index=False)
print(res_df.sort_values("pval")[["metric","variable","stat","pval","fdr"]].to_string(index=False))

# %% ─── Patient-level aggregation and Cox PH ─────────────────────────────────

print("\nPatient-level Cox PH (univariate + Gleason-adjusted)...")
roi_df["pat_id"] = clinical.loc[roi_df.index.intersection(clinical.index), "pat_id"]
pat_shannon = roi_df.groupby("pat_id")[metric_cols].mean()

pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[["psa_progr", "psa_progr_time", "clinical_progr",
                         "clinical_progr_time", "gleason_grp"]]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")
pat_data = pat_shannon.join(pat_clin, how="inner")

cox_results: list[dict] = []
for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    sub = pat_data[[duration_col, event_col, "gleason_grp"] + metric_cols].dropna(
        subset=[duration_col, event_col, "gleason_grp"]
    )
    sub = sub[sub[duration_col] > 0]
    print(f"  {endpoint}: n={len(sub)}, events={int(sub[event_col].sum())}")

    for col in metric_cols:
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
                "n": len(cph_df), "n_events": int(cph_df["E"].sum()),
            })
        except Exception:
            pass

cox_df = pd.DataFrame(cox_results)
if len(cox_df):
    for ep in cox_df["endpoint"].unique():
        mask = cox_df["endpoint"] == ep
        _, fdr, _, _ = multipletests(cox_df.loc[mask, "p_adj"].fillna(1), method="fdr_bh")
        cox_df.loc[mask, "fdr_adj"] = fdr
    cox_df.to_csv(OUT_TAB / "restricted_shannon_cox.csv", index=False)
    print("\nCox results:")
    print(cox_df.sort_values("p_adj")[
        ["metric", "endpoint", "hr_uni", "p_uni", "hr_adj", "ci_low_adj", "ci_high_adj", "p_adj"]
    ].to_string(index=False))

# %% ─── Figure: global vs non-immune Shannon scatter ─────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, hue_col, title in zip(axes,
    ["inflammation", "gleason_grp"],
    ["Inflammation", "Gleason GG"]):
    sub = roi_clin.copy()
    sub["h_global"] = roi_df["h_global"]
    sub["h_non_immune"] = roi_df["h_non_immune"]
    sub = sub.dropna(subset=["h_global", "h_non_immune", hue_col])
    if pd.api.types.is_numeric_dtype(sub[hue_col]) or sub[hue_col].dtype == float:
        sc = ax.scatter(sub["h_global"], sub["h_non_immune"],
                        c=pd.to_numeric(sub[hue_col], errors="coerce"), cmap="RdYlBu_r",
                        s=10, alpha=0.6)
        plt.colorbar(sc, ax=ax, label=title)
    else:
        for cat in sorted(sub[hue_col].unique()):
            m = sub[sub[hue_col] == cat]
            ax.scatter(m["h_global"], m["h_non_immune"], label=str(cat), s=10, alpha=0.6)
        ax.legend(fontsize=7, title=title)
    ax.set_xlabel("Global Shannon (all cells)")
    ax.set_ylabel("Non-immune Shannon")
    ax.set_title(f"Shannon decomposition — {title}")
plt.tight_layout()
fig.savefig(OUT_FIG / "global_vs_nonimmune_shannon.pdf", bbox_inches="tight")
plt.close()

print("\nRQ17 complete. Outputs →", OUT_TAB, "and", OUT_FIG)
