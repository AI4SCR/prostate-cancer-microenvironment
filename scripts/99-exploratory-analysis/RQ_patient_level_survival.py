"""
Patient-level spatial score aggregation and survival analysis.

Aggregates ROI-level spatial features to patient level (max-pool per draft paper approach)
and runs Kaplan-Meier and univariate Cox PH analyses for key spatial scores.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/survival")
OUT_TAB = Path("output/tables/survival")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load clinical (patient-level) ────────────────────────────────────────

print("Loading data...")
clinical = pd.read_parquet(BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet")

# Patient-level clinical: one row per patient
# All ROIs of the same patient share the same clinical outcome
ENDPOINT_COLS = ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time",
                  "disease_progr", "disease_progr_time"]
STRAT_COLS = ["gleason_grp", "stromogenic_smc_loss_reactive_stroma_present", "inflammation",
              "cribriform", "age_at_surgery", "psa_at_surgery", "d_amico_risk"]

pat_clinical = clinical.groupby("pat_id").first()[ENDPOINT_COLS + STRAT_COLS].copy()
print(f"  {len(pat_clinical)} patients")
print(f"  PSA progression events: {pat_clinical.psa_progr.sum()}")
print(f"  Clinical progression events: {pat_clinical.clinical_progr.sum()}")

# %% ─── Load ROI-level spatial scores ─────────────────────────────────────────

# TLS scores (from RQ3)
tls_raw = pd.read_csv("output/tables/RQ3/roi_tls_scores.csv", index_col=0)
tls_scores = tls_raw[["n_bcells", "bcell_proportion", "tls_score", "tls_thelper_score",
                        "tls_treg_score"]].copy()
# Avoid float precision loss: read pat_id from clinical directly using sample_id index
tls_scores["pat_id"] = clinical.loc[tls_scores.index.intersection(clinical.index), "pat_id"]

# CAF proximity scores (from RQ4)
caf_raw = pd.read_csv("output/tables/RQ4/roi_caf_proximity_profiles.csv", index_col=0)
caf_scores = caf_raw[["caf1_plus_periglandular_fraction", "caf2_ar_periglandular_fraction",
                        "enrichment_caf1plus_luminal", "enrichment_caf2ar_luminal",
                        "prop_caf1_plus", "prop_caf2_ar", "prop_luminal"]].copy()

# Cell composition scores (from RQ1 — recompute directly here for reliability)
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
tumor_ids = set(clinical.index[clinical.is_tumor == "yes"])
cell_meta = pd.concat(
    [pd.read_parquet(f, engine="fastparquet") for f in meta_files], ignore_index=False
)
cell_meta = cell_meta.loc[cell_meta.index.get_level_values("sample_id").isin(tumor_ids)]
meta_counts = cell_meta.groupby(["sample_id", "meta_label"]).size().unstack(fill_value=0)
meta_props = meta_counts.div(meta_counts.sum(axis=1), axis=0)

# Add pat_id to meta_props
meta_props_with_pat = meta_props.join(clinical[["pat_id"]])
meta_props_with_pat["pat_id"] = meta_props_with_pat["pat_id"].astype(str)

# %% ─── Combine ROI-level scores ──────────────────────────────────────────────

roi_scores = tls_scores.join(caf_scores, how="outer")
roi_scores = roi_scores.join(meta_props_with_pat, how="left", rsuffix="_comp")
# use pat_id from tls_scores (both should be same)
roi_scores["pat_id"] = roi_scores["pat_id"].fillna(roi_scores.get("pat_id_comp", pd.Series(dtype=str)))
assert "pat_id" in roi_scores.columns

# %% ─── Patient-level max-pooling ────────────────────────────────────────────

SCORE_COLS = [
    "bcell_proportion", "tls_score", "tls_thelper_score", "tls_treg_score",
    "caf1_plus_periglandular_fraction", "caf2_ar_periglandular_fraction",
    "enrichment_caf1plus_luminal", "enrichment_caf2ar_luminal",
    "prop_caf1_plus", "prop_caf2_ar",
    "stromal-CAF2-AR+", "stromal-CAF1-CD105+", "immune-B-cells", "immune-T-cells", "immune-myeloid",
]
# Keep only score columns + pat_id
cols_to_pool = [c for c in SCORE_COLS if c in roi_scores.columns]

roi_for_pool = roi_scores[cols_to_pool + ["pat_id"]].dropna(subset=["pat_id"])
pat_scores = roi_for_pool.groupby("pat_id")[cols_to_pool].max()
print(f"\nPatient-level scores shape: {pat_scores.shape}")

# Join with clinical
pat_all = pat_scores.join(pat_clinical, how="inner")
print(f"Patients with scores + clinical: {len(pat_all)}")
print(f"PSA events: {pat_all.psa_progr.sum()}")
print(f"Clinical progression events: {pat_all.clinical_progr.sum()}")

pat_all.to_csv(OUT_TAB / "patient_level_scores.csv")

# %% ─── Correlations with Gleason ────────────────────────────────────────────

print("\nCorrelations with max Gleason grade group per patient:")
for col in cols_to_pool:
    sub = pat_all.dropna(subset=[col, "gleason_grp"])
    if len(sub) < 20:
        continue
    rho, pval = spearmanr(sub["gleason_grp"].values, sub[col].values)
    if pval < 0.1:
        print(f"  {col}: rho={rho:.3f}, p={pval:.4f} *")

# %% ─── Kaplan-Meier stratified by spatial scores ─────────────────────────────

print("\nKaplan-Meier analyses...")

survival_endpoints = [
    ("psa_progr", "psa_progr_time", "PSA recurrence"),
    ("clinical_progr", "clinical_progr_time", "Clinical progression"),
]

score_labels = {
    "caf1_plus_periglandular_fraction": "CAF1-CD105+ periglandular fraction",
    "enrichment_caf1plus_luminal": "CAF1+ luminal contact enrichment",
    "bcell_proportion": "B-cell proportion",
    "stromal-CAF2-AR+": "CAF2-AR+ proportion",
    "tls_score": "TLS score",
    "immune-B-cells": "B-cell proportion (composition)",
}

km_results = []
fig, axes = plt.subplots(len(score_labels), len(survival_endpoints),
                          figsize=(12, 5 * len(score_labels)))

for row_idx, (score_col, score_label) in enumerate(score_labels.items()):
    if score_col not in pat_all.columns:
        continue
    for col_idx, (event_col, time_col, endpoint_label) in enumerate(survival_endpoints):
        ax = axes[row_idx, col_idx]

        sub = pat_all.dropna(subset=[score_col, event_col, time_col]).copy()
        if len(sub) < 20 or sub[event_col].sum() < 5:
            ax.set_title(f"Insufficient events (n={sub[event_col].sum()})")
            ax.axis("off")
            continue

        median = sub[score_col].median()
        sub["group"] = (sub[score_col] >= median).map({True: "High", False: "Low"})

        g_low = sub[sub.group == "Low"]
        g_high = sub[sub.group == "High"]

        lr = logrank_test(
            g_low[time_col], g_high[time_col],
            event_observed_A=g_low[event_col],
            event_observed_B=g_high[event_col],
        )

        KaplanMeierFitter().fit(g_low[time_col], g_low[event_col],
                                 label=f"Low (n={len(g_low)})").plot_survival_function(
            ax=ax, ci_show=False, color="steelblue")
        KaplanMeierFitter().fit(g_high[time_col], g_high[event_col],
                                 label=f"High (n={len(g_high)})").plot_survival_function(
            ax=ax, ci_show=False, color="salmon")

        ax.set_title(f"{score_label}\n{endpoint_label} — p={lr.p_value:.3f}", fontsize=8)
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Survival")
        ax.legend(fontsize=7)

        km_results.append({
            "score": score_col, "endpoint": event_col, "n": len(sub),
            "n_events": sub[event_col].sum(), "logrank_p": lr.p_value,
        })
        print(f"  {score_col} vs {event_col}: p={lr.p_value:.4f} (n={len(sub)}, events={sub[event_col].sum()})")

plt.tight_layout()
fig.savefig(OUT_FIG / "km_spatial_scores.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

km_df = pd.DataFrame(km_results).sort_values("logrank_p")
km_df.to_csv(OUT_TAB / "km_logrank_results.csv", index=False)

# %% ─── Univariate Cox PH ─────────────────────────────────────────────────────

print("\nUnivariate Cox PH models...")

cox_results = []
for score_col, score_label in score_labels.items():
    if score_col not in pat_all.columns:
        continue
    for event_col, time_col, endpoint_label in survival_endpoints:
        sub = pat_all.dropna(subset=[score_col, event_col, time_col]).copy()
        if len(sub) < 20 or sub[event_col].sum() < 5:
            continue
        try:
            cph = CoxPHFitter()
            cox_df = sub[[score_col, event_col, time_col]].copy()
            cox_df[score_col] = (cox_df[score_col] - cox_df[score_col].mean()) / (cox_df[score_col].std() + 1e-10)
            cph.fit(cox_df, duration_col=time_col, event_col=event_col, show_progress=False)
            row = cph.summary.loc[score_col]
            cox_results.append({
                "score": score_label, "endpoint": endpoint_label,
                "hazard_ratio": float(np.exp(row["coef"])),
                "hr_lower": float(np.exp(row["coef lower 95%"])),
                "hr_upper": float(np.exp(row["coef upper 95%"])),
                "p": float(row["p"]),
                "n": len(sub), "n_events": int(sub[event_col].sum()),
            })
        except Exception as e:
            print(f"  Cox failed: {score_col}/{endpoint_label}: {e}")

cox_df = pd.DataFrame(cox_results).sort_values("p")
cox_df.to_csv(OUT_TAB / "cox_results.csv", index=False)
print("\nUnivariate Cox results:")
print(cox_df.to_string(index=False))

# %% ─── Forest plot ───────────────────────────────────────────────────────────

if len(cox_df) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, endpoint_label in zip(axes, ["PSA recurrence", "Clinical progression"]):
        sub_cox = cox_df[cox_df.endpoint == endpoint_label].copy()
        if len(sub_cox) == 0:
            ax.axis("off")
            continue
        sub_cox = sub_cox.sort_values("hazard_ratio")
        y_pos = list(range(len(sub_cox)))
        ax.errorbar(
            sub_cox.hazard_ratio, y_pos,
            xerr=[sub_cox.hazard_ratio - sub_cox.hr_lower,
                  sub_cox.hr_upper - sub_cox.hazard_ratio],
            fmt="o", color="black", capsize=4,
        )
        ax.axvline(1.0, color="red", linestyle="--", linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(
            [f"{r['score']}\np={r['p']:.3f}" for _, r in sub_cox.iterrows()],
            fontsize=8,
        )
        ax.set_xlabel("Hazard ratio (95% CI)")
        ax.set_title(f"Univariate Cox PH — {endpoint_label}\n(standardized scores)")

    plt.tight_layout()
    fig.savefig(OUT_FIG / "forest_plot_cox.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nForest plot saved.")

print(f"\nAll outputs saved to {OUT_TAB}/ and {OUT_FIG}/")
