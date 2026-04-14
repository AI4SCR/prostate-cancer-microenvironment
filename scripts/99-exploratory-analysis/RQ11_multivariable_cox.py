"""
RQ11: Multivariable Cox Models — TME Scores Adjusted for Gleason Grade

Tests whether spatial TME scores predict PSA recurrence and clinical progression
independently of Gleason grade.

Candidate scores (from prior analyses):
- global_shannon: mean per-ROI global Shannon diversity (RQ5)
- local_shannon_mean: mean local Shannon diversity (RQ5)
- caf1_periglandular_frac: CAF1-CD105+ fraction within 32px of luminal epithelium (RQ4)
- caf2_ar_plus_prop: CAF2-AR+ cell proportion (from RQ1 / clinical or re-computed)
- bcell_max_cluster: max B-cell cluster size (RQ10)
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ11")
OUT_TAB = Path("output/tables/RQ11")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load clinical data ────────────────────────────────────────────────────

print("Loading clinical data...")
clinical = pd.read_parquet(
    BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet"
)
tumor_ids = set(clinical.index[clinical.is_tumor == "yes"])

pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[
        ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time", "gleason_grp"]
    ]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")
print(f"  {len(pat_clin)} patients")

# %% ─── Load TME scores from prior analyses ────────────────────────────────────

scores: dict[str, pd.Series] = {}

# ── RQ5: Shannon diversity ──
rq5_path = Path("output/tables/RQ5/roi_shannon_diversity.csv")
if rq5_path.exists():
    rq5 = pd.read_csv(rq5_path, index_col=0)
    rq5.index = rq5.index.astype(str)
    # attach pat_id
    rq5["pat_id"] = clinical.loc[rq5.index.intersection(clinical.index), "pat_id"]
    for col in ["h_global", "h_local_mean"]:
        if col in rq5.columns:
            pat_score = rq5.groupby("pat_id")[col].mean()
            scores[col] = pat_score
    print(f"  RQ5: loaded {list(scores.keys())}")
else:
    print(f"  WARNING: {rq5_path} not found — skipping Shannon diversity scores")

# ── RQ4: CAF1 periglandular fraction ──
rq4_path = Path("output/tables/RQ4/roi_caf_proximity_profiles.csv")
if rq4_path.exists():
    rq4 = pd.read_csv(rq4_path, index_col=0)
    rq4.index = rq4.index.astype(str)
    rq4["pat_id"] = clinical.loc[rq4.index.intersection(clinical.index), "pat_id"]
    col = "caf1_plus_periglandular_fraction"
    if col in rq4.columns:
        scores["caf1_periglandular_frac"] = rq4.groupby("pat_id")[col].mean()
        print(f"  RQ4: loaded caf1_periglandular_frac")
    else:
        print(f"  WARNING: {col} not found in RQ4 table (available: {list(rq4.columns)})")
else:
    print(f"  WARNING: {rq4_path} not found — skipping CAF1 periglandular fraction")

# ── RQ10: B-cell max cluster size ──
rq10_path = Path("output/tables/RQ10/roi_bcell_cluster_metrics.csv")
if rq10_path.exists():
    rq10 = pd.read_csv(rq10_path, index_col=0)
    rq10.index = rq10.index.astype(str)
    rq10["pat_id"] = clinical.loc[rq10.index.intersection(clinical.index), "pat_id"]
    scores["bcell_max_cluster"] = rq10.groupby("pat_id")["max_cluster_size"].max()
    print(f"  RQ10: loaded bcell_max_cluster")
else:
    print(f"  WARNING: {rq10_path} not found — skipping B-cell cluster size")

# ── CAF2-AR+ proportion: recompute from per-ROI meta-label data ──
print("  Computing CAF2-AR+ proportion from cell data...")
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
meta_by_id = {f.stem: f for f in meta_files}
sample_ids = sorted(set(meta_by_id) & tumor_ids)

caf2_props: list[dict] = []
for sid in sample_ids:
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    if "sample_id" in meta_df.index.names:
        meta_df = meta_df.xs(sid, level="sample_id") if sid in meta_df.index.get_level_values("sample_id") else meta_df
    lbl = meta_df["meta_label"].values
    n = len(lbl)
    if n == 0:
        continue
    prop = (lbl == "stromal-CAF2-AR+").sum() / n
    caf2_props.append({"sample_id": sid, "caf2_ar_plus_prop": prop})

caf2_df = pd.DataFrame(caf2_props).set_index("sample_id")
caf2_df["pat_id"] = clinical.loc[caf2_df.index.intersection(clinical.index), "pat_id"]
scores["caf2_ar_plus_prop"] = caf2_df.groupby("pat_id")["caf2_ar_plus_prop"].mean()
print(f"  CAF2-AR+ proportion: {len(scores['caf2_ar_plus_prop'])} patients")

# %% ─── Assemble patient-level data ───────────────────────────────────────────

score_df = pd.DataFrame(scores)
print(f"\nScore matrix: {score_df.shape[0]} patients × {score_df.shape[1]} scores")
print(f"  Scores: {list(score_df.columns)}")

pat_data = score_df.join(pat_clin, how="inner")
print(f"  After join with clinical: {len(pat_data)} patients")

# Standardize scores to z-scores for comparability
score_cols = list(scores.keys())
for col in score_cols:
    mu, sd = pat_data[col].mean(), pat_data[col].std()
    if sd > 0:
        pat_data[col + "_z"] = (pat_data[col] - mu) / sd
    else:
        pat_data[col + "_z"] = 0.0

z_cols = [c + "_z" for c in score_cols]

# %% ─── Multivariable Cox: score + Gleason ────────────────────────────────────

print("\nRunning multivariable Cox (score + Gleason) per endpoint...")
results: list[dict] = []

for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    sub = pat_data[[duration_col, event_col, "gleason_grp"] + z_cols].dropna(subset=[duration_col, event_col, "gleason_grp"])
    sub = sub[sub[duration_col] > 0]
    n_events = int(sub[event_col].sum())
    print(f"  {endpoint}: n={len(sub)}, events={n_events}")

    for zcol in z_cols:
        orig_col = zcol.replace("_z", "")
        cph_df = sub[[duration_col, event_col, "gleason_grp", zcol]].dropna()
        cph_df = cph_df.rename(columns={duration_col: "T", event_col: "E"})
        if cph_df[zcol].std() < 1e-8:
            continue
        try:
            cph = CoxPHFitter()
            cph.fit(cph_df, duration_col="T", event_col="E", show_progress=False)
            s = cph.summary.loc[zcol]
            g = cph.summary.loc["gleason_grp"]
            results.append({
                "score": orig_col,
                "endpoint": endpoint,
                "n": len(cph_df),
                "n_events": n_events,
                "hr_score": np.exp(s["coef"]),
                "ci_low_score": np.exp(s["coef lower 95%"]),
                "ci_high_score": np.exp(s["coef upper 95%"]),
                "p_score": s["p"],
                "hr_gleason": np.exp(g["coef"]),
                "p_gleason": g["p"],
            })
        except Exception as e:
            print(f"    Failed for {zcol}/{endpoint}: {e}")

res_df = pd.DataFrame(results)
res_df.to_csv(OUT_TAB / "multivariable_cox_results.csv", index=False)
print("\nMultivariable Cox results (score coefficient, adjusted for Gleason):")
print(res_df[["score", "endpoint", "hr_score", "ci_low_score", "ci_high_score", "p_score", "hr_gleason", "p_gleason"]].to_string(index=False))

# %% ─── Full model: Gleason + all scores ─────────────────────────────────────

print("\nFull model: Gleason + all scores simultaneously...")
full_results: list[dict] = []

for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    cols_needed = [duration_col, event_col, "gleason_grp"] + z_cols
    sub = pat_data[cols_needed].dropna()
    sub = sub[sub[duration_col] > 0]
    sub = sub.rename(columns={duration_col: "T", event_col: "E"})
    n_events = int(sub["E"].sum())
    # rule of thumb: need ≥10 events per covariate; skip if underpowered
    n_cov = len(z_cols) + 1
    if n_events < n_cov * 10:
        print(f"  {endpoint}: n_events={n_events}, {n_cov} covariates — underpowered, skipping")
        continue
    try:
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(sub, duration_col="T", event_col="E", show_progress=False)
        print(f"\n  Full model — {endpoint} (n={len(sub)}, events={n_events}):")
        for coef in cph.summary.index:
            s = cph.summary.loc[coef]
            print(f"    {coef:35s}  HR={np.exp(s['coef']):.3f}  p={s['p']:.3f}")
            full_results.append({
                "covariate": coef,
                "endpoint": endpoint,
                "hr": np.exp(s["coef"]),
                "ci_low": np.exp(s["coef lower 95%"]),
                "ci_high": np.exp(s["coef upper 95%"]),
                "p": s["p"],
            })
    except Exception as e:
        print(f"  Full model failed for {endpoint}: {e}")

if full_results:
    pd.DataFrame(full_results).to_csv(OUT_TAB / "full_model_cox_results.csv", index=False)

# %% ─── Forest plot: adjusted HRs ─────────────────────────────────────────────

for ep in res_df["endpoint"].unique():
    ep_df = res_df[res_df["endpoint"] == ep].copy().reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    y = range(len(ep_df))
    colors = ["tomato" if hr > 1 else "steelblue" for hr in ep_df["hr_score"]]
    ax.scatter(ep_df["hr_score"], list(y), color=colors, zorder=3, s=50)
    ax.hlines(list(y), ep_df["ci_low_score"], ep_df["ci_high_score"], color="gray", linewidth=1)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(ep_df["score"], fontsize=9)
    ax.set_xlabel("Hazard Ratio per SD (95% CI)\nadjusted for Gleason grade")
    ax.set_title(f"TME scores adjusted for Gleason — {ep}")
    # annotate p-values
    for i, row in ep_df.iterrows():
        ax.text(ep_df["ci_high_score"].max() * 1.02, i, f"p={row['p_score']:.3f}", va="center", fontsize=8)
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"forest_adjusted_{ep}.pdf", bbox_inches="tight")
    plt.close()
    print(f"  Saved forest plot → {OUT_FIG / f'forest_adjusted_{ep}.pdf'}")

print("\nRQ11 complete. Outputs saved to", OUT_TAB, "and", OUT_FIG)
