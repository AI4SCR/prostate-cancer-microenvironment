"""
RQ16: CD8+ T-Cell Spatial Exclusion and T-Cell Subtype Balance

Computes clinically relevant immune metrics at the label level:
  1. CD8:FoxP3 ratio per ROI (effector:regulatory T-cell balance)
  2. CD8+ spatial exclusion: mean distance of cytotoxic T cells to tumor epithelium
  3. Fraction of CD8+ cells in direct contact (radius-32) with tumor/luminal epithelium
  4. Immune phenotype classification: inflamed / excluded / desert

Tests all metrics vs Gleason, stromogenic status, PSA recurrence, and clinical progression.
Gleason-adjusted Cox PH for metrics that survive univariate testing.
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
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ16")
OUT_TAB = Path("output/tables/RQ16")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

RADIUS = 32  # ~32 µm neighborhood

# Fine-grained label names
CD8_LABEL   = "immune-T-cells_cytotoxic(CD3+CD8a+)"
FOXP3_LABEL = "immune-T-cells_regulatory(CD3+CD4+FoxP3+)"
THELPER_LABEL = "immune-T-cells_helper(CD3+CD4+)"
TUMOR_LABELS = {"epithelial-tumor", "epithelial-(ERG+CD44+)",
                "epithelial-luminal(ERG+)", "epithelial-luminal(ERG+p53+)"}
LUMINAL_LABEL = "epithelial-luminal"
ALL_EPI_LABELS = TUMOR_LABELS | {LUMINAL_LABEL, "epithelial-luminal(Ki67+)"}

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

# %% ─── Per-ROI immune metrics ────────────────────────────────────────────────

print("Computing CD8+ spatial exclusion metrics per ROI...")
roi_metrics: list[dict] = []

for sid in sample_ids:
    meta_df = pd.read_parquet(meta_by_id[sid], engine="fastparquet")
    spat_df = pd.read_parquet(spatial_by_id[sid], engine="fastparquet")

    if "sample_id" in meta_df.index.names:
        meta_df = meta_df.xs(sid, level="sample_id") if sid in meta_df.index.get_level_values("sample_id") else meta_df
    if "sample_id" in spat_df.index.names:
        spat_df = spat_df.xs(sid, level="sample_id") if sid in spat_df.index.get_level_values("sample_id") else spat_df

    idx = meta_df.index.intersection(spat_df.index)
    labels  = meta_df.loc[idx, "label"].to_numpy()
    coords  = spat_df.loc[idx, ["x", "y"]].to_numpy(dtype=np.float32)
    n_total = len(labels)

    # Cell counts
    n_cd8   = (labels == CD8_LABEL).sum()
    n_foxp3 = (labels == FOXP3_LABEL).sum()
    n_th    = (labels == THELPER_LABEL).sum()
    epi_mask   = np.isin(labels, list(ALL_EPI_LABELS))
    tumor_mask = np.isin(labels, list(TUMOR_LABELS))
    n_epi   = epi_mask.sum()

    row: dict = {
        "sample_id": sid,
        "n_total": n_total,
        "n_cd8": int(n_cd8),
        "n_foxp3": int(n_foxp3),
        "n_thelper": int(n_th),
        "cd8_prop": n_cd8 / n_total if n_total > 0 else 0.0,
        "foxp3_prop": n_foxp3 / n_total if n_total > 0 else 0.0,
        "cd8_foxp3_ratio": n_cd8 / (n_foxp3 + 1),  # +1 pseudocount
        "cd8_thelper_ratio": n_cd8 / (n_th + 1),
    }

    # Spatial: CD8 distance to nearest epithelial cell
    cd8_coords = coords[labels == CD8_LABEL]
    epi_coords = coords[epi_mask]
    tumor_coords = coords[tumor_mask]

    if len(cd8_coords) >= 2 and len(epi_coords) >= 2:
        tree_epi = cKDTree(epi_coords)
        dists, _ = tree_epi.query(cd8_coords, k=1)
        row["cd8_mean_dist_to_epi"]   = float(dists.mean())
        row["cd8_median_dist_to_epi"] = float(np.median(dists))

        # Fraction of CD8+ cells within radius-32 of any epithelial cell
        row["cd8_frac_epi_contact"] = float((dists <= RADIUS).sum() / len(cd8_coords))
    else:
        row["cd8_mean_dist_to_epi"]   = np.nan
        row["cd8_median_dist_to_epi"] = np.nan
        row["cd8_frac_epi_contact"]   = np.nan

    if len(cd8_coords) >= 2 and len(tumor_coords) >= 2:
        tree_tumor = cKDTree(tumor_coords)
        dists_t, _ = tree_tumor.query(cd8_coords, k=1)
        row["cd8_mean_dist_to_tumor"]   = float(dists_t.mean())
        row["cd8_frac_tumor_contact"]   = float((dists_t <= RADIUS).sum() / len(cd8_coords))
    else:
        row["cd8_mean_dist_to_tumor"]   = np.nan
        row["cd8_frac_tumor_contact"]   = np.nan

    roi_metrics.append(row)

roi_df = pd.DataFrame(roi_metrics).set_index("sample_id")
roi_df.to_csv(OUT_TAB / "roi_cd8_metrics.csv")
print(f"  {len(roi_df)} ROIs computed")
print(f"  ROIs with ≥2 CD8+ cells: {(roi_df['n_cd8'] >= 2).sum()}")
print(f"  Mean CD8:FoxP3 ratio: {roi_df['cd8_foxp3_ratio'].median():.2f} (median)")
print(f"  Mean CD8 fraction in epi contact: {roi_df['cd8_frac_epi_contact'].mean():.3f}")

# %% ─── Immune phenotype classification ───────────────────────────────────────

# Inflamed:  CD8+ cells present AND infiltrating epithelium (frac_epi_contact > 0.2)
# Excluded:  CD8+ cells present but NOT infiltrating (frac_epi_contact ≤ 0.2)
# Desert:    Very few CD8+ cells (cd8_prop < 1st quartile)
cd8_q1 = roi_df["cd8_prop"].quantile(0.25)
epi_contact_threshold = 0.2

def classify_immune_phenotype(row: pd.Series) -> str:
    if row["cd8_prop"] <= cd8_q1:
        return "desert"
    if pd.isna(row["cd8_frac_epi_contact"]):
        return "desert"
    if row["cd8_frac_epi_contact"] > epi_contact_threshold:
        return "inflamed"
    return "excluded"

roi_df["immune_phenotype"] = roi_df.apply(classify_immune_phenotype, axis=1)
print(f"\n  Immune phenotype distribution:\n{roi_df['immune_phenotype'].value_counts()}")

# %% ─── Clinical associations ─────────────────────────────────────────────────

print("\nTesting CD8 metrics vs clinical variables...")
roi_clin = clinical.loc[clinical.index.isin(roi_df.index)].copy()

metric_cols = [
    "cd8_prop", "foxp3_prop", "cd8_foxp3_ratio",
    "cd8_mean_dist_to_epi", "cd8_frac_epi_contact",
    "cd8_mean_dist_to_tumor", "cd8_frac_tumor_contact",
]

results: list[dict] = []
for col in metric_cols:
    vals = roi_df.loc[roi_clin.index, col].dropna()

    # Gleason
    gl = pd.to_numeric(roi_clin.loc[vals.index, "gleason_grp"].replace("", np.nan), errors="coerce").dropna()
    common = vals.index.intersection(gl.index)
    if len(common) > 20:
        rho, p = spearmanr(vals.loc[common], gl.loc[common])
        results.append({"metric": col, "variable": "gleason_spearman", "stat": rho, "pval": p, "n": len(common)})

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
            results.append({"metric": col, "variable": var, "stat": rbc, "pval": p,
                             "n": len(g1) + len(g2)})

# Immune phenotype vs clinical
for var in ["inflammation", "cribriform", "stromogenic"]:
    if var not in roi_clin.columns:
        continue
    ct = pd.crosstab(roi_df.loc[roi_clin.index, "immune_phenotype"], roi_clin[var])
    print(f"\n  Immune phenotype × {var}:\n{ct}")

res_df = pd.DataFrame(results)
_, fdr, _, _ = multipletests(res_df["pval"].fillna(1), method="fdr_bh")
res_df["fdr"] = fdr
res_df.to_csv(OUT_TAB / "roi_cd8_associations.csv", index=False)

sig = res_df[res_df["fdr"] < 0.1].sort_values("fdr")
print(f"\n  {len(sig)} significant ROI-level associations (FDR<0.1):")
print(sig[["metric", "variable", "stat", "pval", "fdr", "n"]].to_string(index=False))

# %% ─── Patient-level aggregation ─────────────────────────────────────────────

print("\nAggregating to patient level...")
roi_df["pat_id"] = clinical.loc[roi_df.index.intersection(clinical.index), "pat_id"]

pat_metrics = roi_df.groupby("pat_id")[metric_cols].mean()
# Also max CD8:FoxP3 (captures hottest ROI per patient)
pat_metrics["cd8_foxp3_ratio_max"] = roi_df.groupby("pat_id")["cd8_foxp3_ratio"].max()

pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[["psa_progr", "psa_progr_time", "clinical_progr",
                         "clinical_progr_time", "gleason_grp"]]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")
pat_data = pat_metrics.join(pat_clin, how="inner")
print(f"  {len(pat_data)} patients")

# %% ─── Cox PH: univariate and Gleason-adjusted ───────────────────────────────

print("\nCox PH (univariate + Gleason-adjusted)...")
all_metric_cols = list(pat_metrics.columns)
cox_results: list[dict] = []

for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    sub = pat_data[[duration_col, event_col, "gleason_grp"] + all_metric_cols].dropna(
        subset=[duration_col, event_col, "gleason_grp"]
    )
    sub = sub[sub[duration_col] > 0]
    print(f"  {endpoint}: n={len(sub)}, events={int(sub[event_col].sum())}")

    for col in all_metric_cols:
        cph_df = sub[[duration_col, event_col, col]].dropna().rename(
            columns={duration_col: "T", event_col: "E", col: "score"}
        )
        if cph_df["score"].std() < 1e-8 or len(cph_df) < 30:
            continue
        try:
            # Univariate
            cph = CoxPHFitter()
            cph.fit(cph_df, duration_col="T", event_col="E", show_progress=False)
            s = cph.summary.loc["score"]
            # Gleason-adjusted
            cph_adj_df = sub[[duration_col, event_col, "gleason_grp", col]].dropna().rename(
                columns={duration_col: "T", event_col: "E", col: "score"}
            )
            cph_adj = CoxPHFitter()
            cph_adj.fit(cph_adj_df, duration_col="T", event_col="E", show_progress=False)
            s_adj = cph_adj.summary.loc["score"]
            cox_results.append({
                "metric": col, "endpoint": endpoint,
                "hr_uni": np.exp(s["coef"]),
                "p_uni": s["p"],
                "hr_adj": np.exp(s_adj["coef"]),
                "ci_low_adj": np.exp(s_adj["coef lower 95%"]),
                "ci_high_adj": np.exp(s_adj["coef upper 95%"]),
                "p_adj": s_adj["p"],
                "n": len(cph_df),
                "n_events": int(cph_df["E"].sum()),
            })
        except Exception as e:
            pass

cox_df = pd.DataFrame(cox_results)
if len(cox_df):
    for ep in cox_df["endpoint"].unique():
        mask = cox_df["endpoint"] == ep
        _, fdr, _, _ = multipletests(cox_df.loc[mask, "p_adj"].fillna(1), method="fdr_bh")
        cox_df.loc[mask, "fdr_adj"] = fdr
    cox_df.to_csv(OUT_TAB / "cd8_cox_results.csv", index=False)
    print("\nCox results (sorted by Gleason-adjusted p):")
    print(cox_df.sort_values("p_adj")[
        ["metric", "endpoint", "hr_uni", "p_uni", "hr_adj", "ci_low_adj", "ci_high_adj", "p_adj", "fdr_adj"]
    ].to_string(index=False))

# %% ─── Figures ───────────────────────────────────────────────────────────────

# 1. CD8 distance to epithelium vs Gleason
gleason_map = {1: "GG1", 2: "GG2", 3: "GG3", 4: "GG4", 5: "GG5"}
roi_clin2 = roi_clin.copy()
roi_clin2["cd8_mean_dist_to_epi"] = roi_df["cd8_mean_dist_to_epi"]
roi_clin2["cd8_frac_epi_contact"] = roi_df["cd8_frac_epi_contact"]
roi_clin2["cd8_foxp3_ratio"] = roi_df["cd8_foxp3_ratio"]
roi_clin2["gl_label"] = pd.to_numeric(roi_clin2["gleason_grp"], errors="coerce").map(gleason_map)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for ax, (col, ylabel) in zip(axes, [
    ("cd8_mean_dist_to_epi", "CD8+ mean dist to epithelium (px)"),
    ("cd8_frac_epi_contact", "Fraction CD8+ in epi contact"),
    ("cd8_foxp3_ratio", "CD8:FoxP3 ratio"),
]):
    sub = roi_clin2.dropna(subset=["gl_label", col[0]])
    sns.boxplot(data=sub, x="gl_label", y=col[0],
                order=["GG1","GG2","GG3","GG4","GG5"], ax=ax, showfliers=False)
    ax.set_xlabel("Gleason Grade Group")
    ax.set_ylabel(col[1])
    ax.set_title(col[1])
plt.tight_layout()
fig.savefig(OUT_FIG / "cd8_metrics_vs_gleason.pdf", bbox_inches="tight")
plt.close()

# 2. Forest plot: Gleason-adjusted Cox
for ep in ["psa_recurrence", "clinical_progression"]:
    ep_df = cox_df[cox_df["endpoint"] == ep].sort_values("hr_adj")
    if ep_df.empty:
        continue
    fig, ax = plt.subplots(figsize=(7, 6))
    y = range(len(ep_df))
    colors = ["tomato" if hr > 1 else "steelblue" for hr in ep_df["hr_adj"]]
    ax.scatter(ep_df["hr_adj"], list(y), color=colors, zorder=3, s=40)
    ax.hlines(list(y), ep_df["ci_low_adj"], ep_df["ci_high_adj"], color="gray", linewidth=1)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(ep_df["metric"], fontsize=7)
    for i, (_, row) in enumerate(ep_df.iterrows()):
        ax.text(ep_df["ci_high_adj"].max() * 1.02, i,
                f"p={row['p_adj']:.3f}", va="center", fontsize=7)
    ax.set_xlabel("Hazard Ratio per unit (Gleason-adjusted)")
    ax.set_title(f"CD8+ immune metrics — {ep}")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"forest_{ep}.pdf", bbox_inches="tight")
    plt.close()

# 3. KM for most significant metric
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
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"km_{ep}_{col[:30]}.pdf", bbox_inches="tight")
    plt.close()

print("\nRQ16 complete. Outputs →", OUT_TAB, "and", OUT_FIG)
