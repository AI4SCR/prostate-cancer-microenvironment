"""
RQ13: Grade-Stratified Inter-Patient Heterogeneity

Tests whether high-grade disease has more or less inter-patient TME variance than low-grade
(convergence vs divergence hypothesis). Uses Levene/Brown-Forsythe test for homogeneity of
variance across Gleason grade groups.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import levene, spearmanr
from statsmodels.stats.multitest import multipletests
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ13")
OUT_TAB = Path("output/tables/RQ13")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load data ─────────────────────────────────────────────────────────────

print("Loading patient-level CLR composition and clinical data...")
clinical = pd.read_parquet(
    BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet"
)
pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[["gleason_grp", "psa_progr", "psa_progr_time",
                         "clinical_progr", "clinical_progr_time"]]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")

clr_df = pd.read_csv("output/tables/RQ7/patient_mean_clr_composition.csv", index_col=0)
_float_to_str: dict[float, str] = {}
for p in pat_clin.index:
    try:
        _float_to_str[float(p)] = p
    except ValueError:
        pass
clr_df.index = clr_df.index.map(lambda x: _float_to_str.get(x, str(x)))
clr_df.index.name = "pat_id"
cell_type_cols = [c for c in clr_df.columns if c.startswith(
    ("endothelial-", "epithelial-", "immune-", "stromal-", "undefined")
)]
clr_df = clr_df[cell_type_cols]

pat_data = clr_df.join(pat_clin, how="inner").dropna(subset=["gleason_grp"])
pat_data["gleason_grp"] = pat_data["gleason_grp"].astype(int)
print(f"  {len(pat_data)} patients with CLR + Gleason")

cell_types = cell_type_cols
grade_groups = sorted(pat_data["gleason_grp"].unique())
print(f"  Grade groups: {grade_groups}")
print(f"  Patients per grade: {pat_data['gleason_grp'].value_counts().sort_index().to_dict()}")

# %% ─── Between-patient SD per grade group ────────────────────────────────────

print("\nComputing between-patient SD per grade group per cell type...")
sd_records: list[dict] = []
for gg in grade_groups:
    sub = pat_data[pat_data["gleason_grp"] == gg][cell_types]
    n = len(sub)
    for ct in cell_types:
        sd_records.append({"gleason_grp": gg, "cell_type": ct, "n": n,
                            "sd": sub[ct].std(), "mean": sub[ct].mean()})

sd_df = pd.DataFrame(sd_records)
sd_df.to_csv(OUT_TAB / "grade_stratified_variance.csv", index=False)

# Coefficient of variation (SD / |mean|) to normalize across cell types
sd_df["cv"] = sd_df["sd"] / (sd_df["mean"].abs() + 1e-8)

# %% ─── Levene's test: variance homogeneity across grade groups ──────────────

print("Running Levene's test for variance homogeneity across grade groups...")
levene_results: list[dict] = []

for ct in cell_types:
    groups = [pat_data.loc[pat_data["gleason_grp"] == gg, ct].dropna().values
              for gg in grade_groups if len(pat_data[pat_data["gleason_grp"] == gg]) >= 5]
    if len(groups) < 2:
        continue
    stat, p = levene(*groups, center="median")  # Brown-Forsythe variant
    # Direction: Spearman rho between grade and SD (is variance increasing or decreasing?)
    grade_sd = sd_df[sd_df["cell_type"] == ct].set_index("gleason_grp")["sd"]
    gl_vals = grade_sd.index.tolist()
    sd_vals = grade_sd.values.tolist()
    rho, rho_p = spearmanr(gl_vals, sd_vals) if len(gl_vals) >= 3 else (np.nan, np.nan)
    levene_results.append({"cell_type": ct, "levene_stat": stat, "levene_p": p,
                            "rho_sd_vs_grade": rho, "rho_p": rho_p})

lev_df = pd.DataFrame(levene_results)
_, fdr, _, _ = multipletests(lev_df["levene_p"].fillna(1), method="fdr_bh")
lev_df["levene_fdr"] = fdr
lev_df.sort_values("levene_p", inplace=True)
lev_df.to_csv(OUT_TAB / "levene_results.csv", index=False)

sig = lev_df[lev_df["levene_fdr"] < 0.1]
print(f"  {len(sig)} cell types with significant variance heterogeneity across grades (FDR<0.1):")
print(lev_df.head(10)[["cell_type", "levene_stat", "levene_p", "levene_fdr",
                         "rho_sd_vs_grade", "rho_p"]].to_string(index=False))

# %% ─── Figure 1: SD per grade group (selected cell types) ───────────────────

# Pick cell types with top Levene statistics
top_cts = lev_df.head(9)["cell_type"].tolist()
pivot_sd = sd_df[sd_df["cell_type"].isin(top_cts)].pivot(
    index="gleason_grp", columns="cell_type", values="sd"
)

fig, axes = plt.subplots(3, 3, figsize=(12, 10))
for ax, ct in zip(axes.ravel(), top_cts):
    sub = sd_df[sd_df["cell_type"] == ct].sort_values("gleason_grp")
    ax.bar(sub["gleason_grp"].astype(str), sub["sd"])
    ax.set_title(ct[:28], fontsize=8)
    ax.set_xlabel("Gleason GG")
    ax.set_ylabel("Between-patient SD")
    rho_row = lev_df[lev_df["cell_type"] == ct].iloc[0]
    ax.set_title(f"{ct[:24]}\nrho={rho_row['rho_sd_vs_grade']:.2f}", fontsize=7)
plt.suptitle("Between-patient SD per Gleason grade group (top 9 by Levene test)", y=1.01)
plt.tight_layout()
fig.savefig(OUT_FIG / "sd_per_grade_top9.pdf", bbox_inches="tight")
plt.close()

# %% ─── Figure 2: Heatmap of SD per grade × cell type ────────────────────────

pivot_full = sd_df.pivot(index="cell_type", columns="gleason_grp", values="sd")
fig, ax = plt.subplots(figsize=(8, 9))
sns.heatmap(pivot_full, cmap="YlOrRd", ax=ax,
            xticklabels=[f"GG{g}" for g in pivot_full.columns],
            yticklabels=[c[:30] for c in pivot_full.index])
ax.set_title("Between-patient CLR SD per Gleason grade group")
ax.set_xlabel("Gleason Grade Group")
ax.set_ylabel("Cell type")
plt.tight_layout()
fig.savefig(OUT_FIG / "sd_heatmap.pdf", bbox_inches="tight")
plt.close()

# %% ─── Figure 3: Overall TME variance (mean SD across all cell types) ────────

mean_sd_per_grade = sd_df.groupby("gleason_grp")["sd"].mean().reset_index()
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(mean_sd_per_grade["gleason_grp"].astype(str), mean_sd_per_grade["sd"])
ax.set_xlabel("Gleason Grade Group")
ax.set_ylabel("Mean between-patient CLR SD (all cell types)")
ax.set_title("Overall TME inter-patient variance by Gleason grade")
plt.tight_layout()
fig.savefig(OUT_FIG / "mean_sd_per_grade.pdf", bbox_inches="tight")
plt.close()

rho_overall, p_overall = spearmanr(mean_sd_per_grade["gleason_grp"], mean_sd_per_grade["sd"])
print(f"\n  Overall mean SD vs grade: Spearman rho={rho_overall:.3f}, p={p_overall:.3f}")
print(f"  Direction: {'CONVERGING (variance decreases with grade)' if rho_overall < 0 else 'DIVERGING (variance increases with grade)'}")

print("\nRQ13 complete. Outputs saved to", OUT_TAB, "and", OUT_FIG)
