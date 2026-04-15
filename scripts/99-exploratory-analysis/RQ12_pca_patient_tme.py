"""
RQ12: PCA of Patient-Level TME Space

Characterizes the dominant axes of inter-patient TME variation using PCA on CLR-transformed
patient-level composition profiles. Tests whether PCA scores predict survival. Spatial scores
from prior analyses are passively projected onto the PCA space to interpret axes.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ12")
OUT_TAB = Path("output/tables/RQ12")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load clinical data ────────────────────────────────────────────────────

print("Loading clinical data...")
clinical = pd.read_parquet(
    BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet"
)
pat_clin = (
    clinical[clinical.is_tumor == "yes"]
    .groupby("pat_id")[
        ["psa_progr", "psa_progr_time", "clinical_progr", "clinical_progr_time",
         "gleason_grp", "stromogenic_smc_loss_reactive_stroma_present", "inflammation"]
    ]
    .first()
)
pat_clin["gleason_grp"] = pd.to_numeric(pat_clin["gleason_grp"], errors="coerce")

# %% ─── Load patient-level CLR composition ────────────────────────────────────

print("Loading patient-level CLR composition from RQ7...")
clr_df = pd.read_csv("output/tables/RQ7/patient_mean_clr_composition.csv", index_col=0)
CELL_TYPE_COLS = [c for c in clr_df.columns if c.startswith(
    ("endothelial-", "epithelial-", "immune-", "stromal-", "undefined")
)]
clr_df = clr_df[CELL_TYPE_COLS]
# Rebuild string pat_ids (float precision loss in CSV)
_float_to_str: dict[float, str] = {}
for p in pat_clin.index:
    try:
        _float_to_str[float(p)] = p
    except ValueError:
        pass
clr_df.index = clr_df.index.map(lambda x: _float_to_str.get(x, str(x)))
clr_df.index.name = "pat_id"
print(f"  {clr_df.shape[0]} patients × {clr_df.shape[1]} cell types")

# %% ─── Load spatial scores ───────────────────────────────────────────────────

spatial_scores: dict[str, pd.Series] = {}

# Shannon diversity (RQ5)
rq5 = pd.read_csv("output/tables/RQ5/roi_shannon_diversity.csv", index_col=0)
rq5.index = rq5.index.astype(str)
rq5["pat_id"] = clinical.loc[rq5.index.intersection(clinical.index), "pat_id"]
spatial_scores["h_global"] = rq5.groupby("pat_id")["h_global"].mean()
spatial_scores["h_local_mean"] = rq5.groupby("pat_id")["h_local_mean"].mean()

# CAF proximity (RQ4)
rq4 = pd.read_csv("output/tables/RQ4/roi_caf_proximity_profiles.csv", index_col=0)
rq4.index = rq4.index.astype(str)
rq4_pat_id = clinical.loc[rq4.index.intersection(clinical.index), "pat_id"]
rq4["pat_id"] = rq4_pat_id
spatial_scores["caf1_periglandular_frac"] = rq4.groupby("pat_id")["caf1_plus_periglandular_fraction"].mean()

# B-cell clusters (RQ10)
rq10 = pd.read_csv("output/tables/RQ10/roi_bcell_cluster_metrics.csv", index_col=0)
rq10.index = rq10.index.astype(str)
rq10["pat_id"] = clinical.loc[rq10.index.intersection(clinical.index), "pat_id"]
spatial_scores["bcell_max_cluster"] = rq10.groupby("pat_id")["max_cluster_size"].max()

scores_df = pd.DataFrame(spatial_scores)

# %% ─── PCA ───────────────────────────────────────────────────────────────────

print("Running PCA on patient-level CLR profiles...")
# z-score (PCA is sensitive to scale; CLR already centers but not scales variance)
scaler = StandardScaler()
X = scaler.fit_transform(clr_df.values)

pca = PCA()
pca.fit(X)
pc_scores = pd.DataFrame(
    pca.transform(X),
    index=clr_df.index,
    columns=[f"PC{i+1}" for i in range(X.shape[1])],
)

# Explained variance
explained = pd.Series(pca.explained_variance_ratio_, index=pc_scores.columns)
print(f"  PC1: {explained['PC1']:.1%}  PC2: {explained['PC2']:.1%}  PC3: {explained['PC3']:.1%}")
print(f"  Top 5 cumulative: {explained.iloc[:5].cumsum().iloc[-1]:.1%}")

explained.to_csv(OUT_TAB / "explained_variance.csv")

# Loadings
loadings = pd.DataFrame(
    pca.components_.T,
    index=clr_df.columns,
    columns=[f"PC{i+1}" for i in range(X.shape[1])],
)
loadings.to_csv(OUT_TAB / "pca_loadings.csv")

pc_scores.to_csv(OUT_TAB / "patient_pc_scores.csv")

# %% ─── Figure 1: Scree plot ──────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(7, 4))
n_show = 10
ax.bar(range(1, n_show + 1), explained.iloc[:n_show].values * 100)
ax2 = ax.twinx()
ax2.plot(range(1, n_show + 1), explained.iloc[:n_show].cumsum().values * 100,
         "o-", color="tomato")
ax.set_xlabel("Principal Component")
ax.set_ylabel("Variance explained (%)")
ax2.set_ylabel("Cumulative variance (%)", color="tomato")
ax.set_title("PCA scree plot — patient-level CLR composition")
plt.tight_layout()
fig.savefig(OUT_FIG / "scree_plot.pdf", bbox_inches="tight")
plt.close()

# %% ─── Figure 2: Loadings heatmap for PC1–PC4 ───────────────────────────────

fig, axes = plt.subplots(1, 4, figsize=(16, 6))
for ax, pc in zip(axes, ["PC1", "PC2", "PC3", "PC4"]):
    top = loadings[pc].abs().nlargest(10).index
    vals = loadings.loc[top, pc].sort_values()
    colors = ["tomato" if v > 0 else "steelblue" for v in vals]
    ax.barh(range(len(vals)), vals.values, color=colors)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels([v[:22] for v in vals.index], fontsize=7)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_title(f"{pc} ({explained[pc]:.1%})")
    ax.set_xlabel("Loading")
plt.suptitle("Top 10 loadings per PC")
plt.tight_layout()
fig.savefig(OUT_FIG / "loadings_top10.pdf", bbox_inches="tight")
plt.close()

# Full loadings heatmap
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(loadings.iloc[:, :6], cmap="RdBu_r", center=0, ax=ax,
            xticklabels=[f"PC{i+1}" for i in range(6)],
            yticklabels=[l[:28] for l in loadings.index])
ax.set_title("PCA loadings — PC1–6")
plt.tight_layout()
fig.savefig(OUT_FIG / "loadings_heatmap.pdf", bbox_inches="tight")
plt.close()

# %% ─── Figure 3: PC1 vs PC2 biplot colored by clinical variables ─────────────

pat_data = pc_scores.join(pat_clin, how="inner").join(scores_df, how="left")
print(f"  {len(pat_data)} patients with clinical data for biplot")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, (color_col, label) in zip(axes, [
    ("gleason_grp", "Gleason Grade Group"),
    ("stromogenic_smc_loss_reactive_stroma_present", "Stromogenic"),
    ("inflammation", "Inflammation"),
]):
    sub = pat_data.dropna(subset=["PC1", "PC2", color_col])
    if sub[color_col].dtype in [float, int] or pd.api.types.is_numeric_dtype(sub[color_col]):
        sc = ax.scatter(sub["PC1"], sub["PC2"], c=sub[color_col], cmap="RdYlBu_r", s=20, alpha=0.7)
        plt.colorbar(sc, ax=ax, label=label)
    else:
        cats = sorted(sub[color_col].dropna().unique())
        palette = dict(zip(cats, sns.color_palette("Set2", len(cats))))
        for cat in cats:
            m = sub[sub[color_col] == cat]
            ax.scatter(m["PC1"], m["PC2"], label=str(cat), s=20, alpha=0.7, color=palette[cat])
        ax.legend(fontsize=7, title=label)
    ax.set_xlabel(f"PC1 ({explained['PC1']:.1%})")
    ax.set_ylabel(f"PC2 ({explained['PC2']:.1%})")
    ax.set_title(f"PC1 vs PC2 — {label}")

plt.tight_layout()
fig.savefig(OUT_FIG / "pc1_pc2_biplot.pdf", bbox_inches="tight")
plt.close()

# %% ─── Correlate PCs with clinical and spatial variables ─────────────────────

print("\nCorrelating PCs with clinical and spatial scores...")
corr_results: list[dict] = []
test_vars = {
    "gleason_grp": "Gleason GG",
    "h_global": "Global Shannon",
    "h_local_mean": "Local Shannon",
    "caf1_periglandular_frac": "CAF1 periglandular",
    "bcell_max_cluster": "B-cell cluster",
}

for pc in ["PC1", "PC2", "PC3", "PC4"]:
    for var, label in test_vars.items():
        sub = pat_data[[pc, var]].dropna()
        if len(sub) < 20:
            continue
        rho, p = spearmanr(sub[pc], sub[var])
        corr_results.append({"PC": pc, "variable": label, "rho": rho, "p": p, "n": len(sub)})

corr_df = pd.DataFrame(corr_results)
corr_df.to_csv(OUT_TAB / "pc_correlation_with_variables.csv", index=False)
print(corr_df.sort_values("p").head(12).to_string(index=False))

# Heatmap of rho
pivot = corr_df.pivot(index="variable", columns="PC", values="rho")
fig, ax = plt.subplots(figsize=(7, 4))
sns.heatmap(pivot, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax)
ax.set_title("Spearman correlation: PC scores vs clinical/spatial variables")
plt.tight_layout()
fig.savefig(OUT_FIG / "pc_variable_correlations.pdf", bbox_inches="tight")
plt.close()

# %% ─── Cox PH: PC1 and PC2 vs survival ───────────────────────────────────────

print("\nCox PH: PC1/PC2 vs survival endpoints...")
cox_results: list[dict] = []

for endpoint, duration_col, event_col in [
    ("psa_recurrence", "psa_progr_time", "psa_progr"),
    ("clinical_progression", "clinical_progr_time", "clinical_progr"),
]:
    sub = pat_data[[duration_col, event_col] + [f"PC{i+1}" for i in range(4)]].dropna(subset=[duration_col, event_col])
    sub = sub[sub[duration_col] > 0]
    n_events = int(sub[event_col].sum())
    print(f"  {endpoint}: n={len(sub)}, events={n_events}")

    for pc in ["PC1", "PC2", "PC3", "PC4"]:
        cph_df = sub[[duration_col, event_col, pc]].dropna().rename(
            columns={duration_col: "T", event_col: "E", pc: "score"}
        )
        if cph_df["score"].std() < 1e-8:
            continue
        try:
            cph = CoxPHFitter()
            cph.fit(cph_df, duration_col="T", event_col="E", show_progress=False)
            s = cph.summary.loc["score"]
            cox_results.append({
                "PC": pc, "endpoint": endpoint,
                "hr": np.exp(s["coef"]),
                "ci_low": np.exp(s["coef lower 95%"]),
                "ci_high": np.exp(s["coef upper 95%"]),
                "p": s["p"],
                "n": len(cph_df),
                "n_events": n_events,
            })
        except Exception as e:
            print(f"    Cox failed {pc}/{endpoint}: {e}")

cox_df = pd.DataFrame(cox_results)
if len(cox_df):
    cox_df.to_csv(OUT_TAB / "pc_cox_results.csv", index=False)
    print("\n  Cox results:")
    print(cox_df[["PC", "endpoint", "hr", "ci_low", "ci_high", "p"]].sort_values("p").to_string(index=False))

# %% ─── KM for most significant PC ────────────────────────────────────────────

if len(cox_df):
    best = cox_df.sort_values("p").iloc[0]
    pc, ep = best["PC"], best["endpoint"]
    duration_col = "psa_progr_time" if ep == "psa_recurrence" else "clinical_progr_time"
    event_col    = "psa_progr"       if ep == "psa_recurrence" else "clinical_progr"

    sub = pat_data[[duration_col, event_col, pc]].dropna()
    sub = sub[sub[duration_col] > 0]
    med = sub[pc].median()
    hi, lo = sub[sub[pc] >= med], sub[sub[pc] < med]
    lr = logrank_test(hi[duration_col], lo[duration_col], hi[event_col], lo[event_col])

    fig, ax = plt.subplots(figsize=(6, 4))
    KaplanMeierFitter().fit(hi[duration_col], hi[event_col], label=f"High {pc} (n={len(hi)})").plot_survival_function(ax=ax)
    KaplanMeierFitter().fit(lo[duration_col], lo[event_col], label=f"Low {pc} (n={len(lo)})").plot_survival_function(ax=ax)
    ax.set_title(f"KM — {ep}\n{pc} high vs low (log-rank p={lr.p_value:.3f}, HR={best['hr']:.2f})")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    plt.tight_layout()
    fig.savefig(OUT_FIG / f"km_{ep}_{pc}.pdf", bbox_inches="tight")
    plt.close()

print("\nRQ12 complete. Outputs saved to", OUT_TAB, "and", OUT_FIG)
