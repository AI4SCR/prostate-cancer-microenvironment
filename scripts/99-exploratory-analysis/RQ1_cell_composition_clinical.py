"""
RQ1: Cell Type Composition Variation Across Clinical and Pathological Groups

Computes ROI-level cell type proportions, applies CLR transformation, and tests
for systematic differences across Gleason grade, stromogenic status, inflammation,
cribriform architecture, PSA recurrence, and clinical progression.
"""

# %% ─── Setup ────────────────────────────────────────────────────────────────

from pathlib import Path
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import kruskal, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from dotenv import load_dotenv
import os

assert load_dotenv(override=True)
BASE_DIR = Path(os.environ["BASE_DIR"])
OUT_FIG = Path("output/figures/RQ1")
OUT_TAB = Path("output/tables/RQ1")
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_TAB.mkdir(parents=True, exist_ok=True)

# %% ─── Load data ─────────────────────────────────────────────────────────────

print("Loading cell-level metadata...")
meta_files = sorted((BASE_DIR / "02_processed/metadata/filtered-annotated").glob("*.parquet"))
assert len(meta_files) > 0, "No metadata files found"

cell_meta = pd.concat(
    [pd.read_parquet(f, engine="fastparquet") for f in meta_files],
    ignore_index=False,
)
assert cell_meta.index.names == ["sample_id", "object_id"]
print(f"  {len(cell_meta):,} cells across {cell_meta.index.get_level_values('sample_id').nunique()} ROIs")

clinical = pd.read_parquet(BASE_DIR / "02_processed/metadata/clinical.parquet", engine="fastparquet")
assert clinical.index.name == "sample_id"
print(f"  Clinical metadata: {clinical.shape}")

# Filter to tumor ROIs only
tumor_ids = clinical.index[clinical.is_tumor == "yes"].tolist()
print(f"  Tumor ROIs: {len(tumor_ids)}")

cell_meta = cell_meta.loc[cell_meta.index.get_level_values("sample_id").isin(tumor_ids)]
clinical = clinical.loc[tumor_ids]

# %% ─── Compute ROI-level cell type proportions ───────────────────────────────

print("Computing cell type proportions...")

# Fine label proportions
label_counts = (
    cell_meta.groupby(["sample_id", "label"]).size().unstack(fill_value=0)
)
label_props = label_counts.div(label_counts.sum(axis=1), axis=0)

# Meta-label proportions
meta_counts = (
    cell_meta.groupby(["sample_id", "meta_label"]).size().unstack(fill_value=0)
)
meta_props = meta_counts.div(meta_counts.sum(axis=1), axis=0)

# Main group proportions
main_counts = (
    cell_meta.groupby(["sample_id", "main_group"]).size().unstack(fill_value=0)
)
main_props = main_counts.div(main_counts.sum(axis=1), axis=0)

print(f"  Label props shape: {label_props.shape}")
print(f"  Meta-label props shape: {meta_props.shape}")

# %% ─── CLR transformation ────────────────────────────────────────────────────

def clr(df: pd.DataFrame, pseudocount: float = 1e-5) -> pd.DataFrame:
    """Centered log-ratio transform for compositional data."""
    X = df.values + pseudocount
    log_X = np.log(X)
    log_gm = log_X.mean(axis=1, keepdims=True)
    return pd.DataFrame(log_X - log_gm, index=df.index, columns=df.columns)

label_clr = clr(label_props)
meta_clr = clr(meta_props)
main_clr = clr(main_props)

# %% ─── Statistical tests across clinical groups ──────────────────────────────

print("Running statistical tests...")

# Clinical variables: {column: (test_type, group_order_or_None)}
CLINICAL_VARS = {
    "gleason_grp": ("kruskal_spearman", None),
    "stromogenic_smc_loss_reactive_stroma_present": ("mannwhitney", ["no", "yes"]),
    "inflammation": ("mannwhitney", ["no", "yes"]),
    "cribriform": ("mannwhitney", ["no", "yes"]),
    "psa_progr": ("mannwhitney", [0, 1]),
    "clinical_progr": ("mannwhitney", [0, 1]),
}


def test_association(prop_df: pd.DataFrame, clin: pd.DataFrame, clin_col: str, test_type: str) -> pd.DataFrame:
    """Test association between each cell type proportion and a clinical variable."""
    common = prop_df.index.intersection(clin.index)
    X = prop_df.loc[common]
    y = clin.loc[common, clin_col]

    # Drop NAs in clinical variable
    mask = y.notna()
    X = X.loc[mask]
    y = y.loc[mask]

    results = []
    for ct in X.columns:
        vals = X[ct].values
        groups = y.values

        if test_type == "mannwhitney":
            unique_groups = sorted(y.unique())
            assert len(unique_groups) == 2, f"Expected 2 groups, got {unique_groups}"
            g0 = vals[groups == unique_groups[0]]
            g1 = vals[groups == unique_groups[1]]
            stat, pval = mannwhitneyu(g0, g1, alternative="two-sided")
            # rank-biserial correlation as effect size
            n0, n1 = len(g0), len(g1)
            rbc = 1 - (2 * stat) / (n0 * n1)
            results.append({"cell_type": ct, "statistic": stat, "pval": pval, "effect_size": rbc,
                            "effect_type": "rank_biserial", "test": "mannwhitney",
                            "mean_g0": g0.mean(), "mean_g1": g1.mean()})

        elif test_type == "kruskal_spearman":
            groups_arr = y.values.astype(float)
            # Kruskal-Wallis
            unique_g = sorted(y.unique())
            group_data = [vals[groups == g] for g in unique_g]
            try:
                stat_kw, pval_kw = kruskal(*group_data)
            except ValueError:
                stat_kw, pval_kw = np.nan, np.nan
            # Spearman correlation with ordinal gleason
            rho, pval_sp = spearmanr(vals, groups_arr)
            results.append({"cell_type": ct, "statistic": stat_kw, "pval": pval_kw,
                            "effect_size": rho, "effect_type": "spearman_rho", "test": "kruskal+spearman",
                            "spearman_pval": pval_sp})

    df = pd.DataFrame(results)
    df["fdr"] = multipletests(df["pval"].fillna(1.0), method="fdr_bh")[1]
    df = df.sort_values("fdr")
    return df


all_results = {}
for clin_col, (test_type, _) in CLINICAL_VARS.items():
    res = test_association(meta_clr, clinical, clin_col, test_type)
    res["clin_var"] = clin_col
    all_results[clin_col] = res
    n_sig = (res["fdr"] < 0.1).sum()
    print(f"  {clin_col}: {n_sig} significant (FDR < 0.1)")

# Also run on main groups
for clin_col, (test_type, _) in CLINICAL_VARS.items():
    res_main = test_association(main_clr, clinical, clin_col, test_type)
    res_main["clin_var"] = clin_col
    all_results[f"{clin_col}_main"] = res_main

# Save all results
combined = pd.concat(all_results.values(), ignore_index=True)
combined.to_csv(OUT_TAB / "association_results.csv", index=False)
print(f"  Saved {len(combined)} test results")

# %% ─── Figure 1: Heatmap of mean CLR proportion per clinical group ───────────

print("Generating heatmap figure...")

# Use gleason_grp as the primary stratification
gleason_map = {1.0: "GG1", 2.0: "GG2", 3.0: "GG3", 4.0: "GG4", 5.0: "GG5"}
common = meta_clr.index.intersection(clinical.index)
gg = clinical.loc[common, "gleason_grp"].dropna()
gg_ids = gg.index

mean_by_gg = meta_clr.loc[gg_ids].groupby(gg).mean()
mean_by_gg.index = [gleason_map.get(g, str(g)) for g in mean_by_gg.index]

# Row-normalize for display (z-score across Gleason groups)
mean_by_gg_z = (mean_by_gg - mean_by_gg.mean()) / (mean_by_gg.std() + 1e-10)

fig, ax = plt.subplots(figsize=(14, 6))
sns.heatmap(
    mean_by_gg_z.T,
    cmap="RdBu_r",
    center=0,
    linewidths=0.5,
    ax=ax,
    cbar_kws={"label": "z-score (CLR proportion)"},
)
ax.set_title("Mean CLR cell type proportion by Gleason grade group (z-scored across groups)")
ax.set_xlabel("Gleason grade group")
ax.set_ylabel("Cell type (meta-label)")
plt.tight_layout()
fig.savefig(OUT_FIG / "heatmap_clr_by_gleason.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 2: Heatmap across all clinical variables ──────────────────────

binary_vars = [
    ("stromogenic_smc_loss_reactive_stroma_present", "no", "yes", "Stromogenic"),
    ("inflammation", "no", "yes", "Inflammation"),
    ("cribriform", "no", "yes", "Cribriform"),
    ("psa_progr", 0, 1, "PSA recurrence"),
    ("clinical_progr", 0, 1, "Clinical progression"),
]

# Effect size heatmap: rows = cell types, columns = clinical variables
es_matrix = {}
for col, g0, g1, label in binary_vars:
    res = all_results[col]
    es_dict = res.set_index("cell_type")["effect_size"].to_dict()
    es_matrix[label] = es_dict

es_df = pd.DataFrame(es_matrix)
# Add Gleason Spearman rho
gleason_res = all_results["gleason_grp"].set_index("cell_type")["effect_size"]
es_df["Gleason (rho)"] = gleason_res

# Sort rows by mean absolute effect
es_df = es_df.loc[es_df.abs().mean(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(9, 8))
sns.heatmap(
    es_df,
    cmap="RdBu_r",
    center=0,
    vmin=-0.5,
    vmax=0.5,
    linewidths=0.5,
    ax=ax,
    cbar_kws={"label": "Effect size (rank-biserial / Spearman rho)"},
    annot=True,
    fmt=".2f",
    annot_kws={"size": 8},
)
ax.set_title("Cell type composition association with clinical variables")
ax.set_ylabel("Cell type (meta-label)")
plt.tight_layout()
fig.savefig(OUT_FIG / "heatmap_effect_size_all_vars.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 3: Volcano plots for key clinical variables ───────────────────

print("Generating volcano plots...")

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.ravel()

plot_vars = [
    ("gleason_grp", "Gleason grade (Spearman rho)", "spearman"),
    ("stromogenic_smc_loss_reactive_stroma_present", "Stromogenic vs non-stromogenic", "rbc"),
    ("inflammation", "Inflammation vs no inflammation", "rbc"),
    ("cribriform", "Cribriform vs no cribriform", "rbc"),
    ("psa_progr", "PSA recurrence (1 vs 0)", "rbc"),
    ("clinical_progr", "Clinical progression (1 vs 0)", "rbc"),
]

for ax, (clin_col, title, es_col) in zip(axes, plot_vars):
    res = all_results[clin_col]
    x = res["effect_size"].values
    y = -np.log10(res["fdr"].values + 1e-10)
    labels = res["cell_type"].values

    colors = np.where(res["fdr"].values < 0.1, "red", "gray")
    ax.scatter(x, y, c=colors, s=30, alpha=0.8)

    # Label significant points
    sig = res[res["fdr"] < 0.1]
    for _, row in sig.iterrows():
        ax.annotate(
            row["cell_type"].replace("stromal-", "s-").replace("immune-", "i-").replace("epithelial-", "e-"),
            (row["effect_size"], -np.log10(row["fdr"] + 1e-10)),
            fontsize=7,
            xytext=(5, 2),
            textcoords="offset points",
        )

    ax.axhline(-np.log10(0.1), color="red", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_xlabel("Effect size")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(title, fontsize=9)

plt.suptitle("Cell type composition associations — volcano plots", fontsize=12, y=1.01)
plt.tight_layout()
fig.savefig(OUT_FIG / "volcano_all_vars.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Figure 4: Top significant cell types — boxplots ──────────────────────

print("Generating top significant cell type boxplots...")

# Focus on PSA recurrence and Gleason as key clinical variables
for clin_col, label, groups in [
    ("psa_progr", "PSA recurrence", [0, 1]),
    ("gleason_grp", "Gleason grade group", [1.0, 2.0, 3.0, 4.0, 5.0]),
    ("stromogenic_smc_loss_reactive_stroma_present", "Stromogenic", ["no", "yes"]),
]:
    res = all_results[clin_col]
    top_ct = res[res["fdr"] < 0.2].sort_values("fdr").head(6)["cell_type"].tolist()

    if len(top_ct) == 0:
        print(f"  No significant cell types for {clin_col}")
        continue

    common = meta_props.index.intersection(clinical.index)
    plot_df = meta_props.loc[common, top_ct].copy()
    plot_df["clin_var"] = clinical.loc[common, clin_col].values

    plot_df = plot_df.dropna(subset=["clin_var"])
    plot_df["clin_var"] = plot_df["clin_var"].astype(str)

    n_ct = len(top_ct)
    ncols = min(3, n_ct)
    nrows = (n_ct + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    if n_ct == 1:
        axes = [axes]
    else:
        axes = np.array(axes).ravel()

    for ax, ct in zip(axes, top_ct):
        fdr_val = res[res["cell_type"] == ct]["fdr"].values[0]
        sns.boxplot(data=plot_df, x="clin_var", y=ct, ax=ax, palette="Set2", order=[str(g) for g in groups])
        ax.set_title(f"{ct}\nFDR={fdr_val:.3f}", fontsize=8)
        ax.set_xlabel(label)
        ax.set_ylabel("Proportion")

    # Hide empty axes
    for ax in axes[n_ct:]:
        ax.set_visible(False)

    plt.suptitle(f"Top cell types by {label}", fontsize=11)
    plt.tight_layout()
    safe_col = clin_col.replace("/", "_")
    fig.savefig(OUT_FIG / f"boxplots_top_{safe_col}.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved boxplots for {clin_col} ({len(top_ct)} cell types)")

# %% ─── Figure 5: Main group proportions overview ─────────────────────────────

print("Generating main group composition figure...")

common = main_props.index.intersection(clinical.index)
gg = clinical.loc[common, "gleason_grp"].dropna()

main_mean_gg = main_props.loc[gg.index].groupby(gg).mean()
main_mean_gg.index = [gleason_map.get(g, str(g)) for g in main_mean_gg.index]

fig, ax = plt.subplots(figsize=(8, 5))
main_mean_gg.plot(kind="bar", stacked=True, ax=ax, colormap="tab10")
ax.set_title("Mean cell type main group proportions by Gleason grade group")
ax.set_xlabel("Gleason grade group")
ax.set_ylabel("Mean proportion")
ax.legend(title="Main group", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)
plt.tight_layout()
fig.savefig(OUT_FIG / "barplot_main_groups_by_gleason.pdf", dpi=150, bbox_inches="tight")
plt.close(fig)

# %% ─── Summary table ─────────────────────────────────────────────────────────

print("\nGenerating summary table...")

summary_rows = []
for clin_col, res in all_results.items():
    sig = res[res["fdr"] < 0.1]
    if len(sig) == 0:
        continue
    for _, row in sig.iterrows():
        summary_rows.append({
            "clinical_variable": clin_col,
            "cell_type": row["cell_type"],
            "fdr": row["fdr"],
            "effect_size": row["effect_size"],
            "effect_type": row.get("effect_type", ""),
            "pval": row["pval"],
        })

summary_df = pd.DataFrame(summary_rows).sort_values(["clinical_variable", "fdr"])
summary_df.to_csv(OUT_TAB / "significant_associations.csv", index=False)
print(f"  {len(summary_df)} significant associations (FDR < 0.1)")
print("\nTop associations:")
print(summary_df.head(20).to_string(index=False))

print("\nDone. All outputs written to output/figures/RQ1/ and output/tables/RQ1/")
