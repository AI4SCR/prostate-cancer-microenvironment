# Exploratory Analysis Findings Summary

**Dataset:** 190 patients, 476 tumor ROIs, 2.19M cells, 34 cell types / 17 meta-labels, 40 markers (IMC)  
**Date:** 2026-04-13  
**Scripts:** `scripts/99-exploratory-analysis/RQ1–4*.py` and `RQ_patient_level_survival.py`  
**Outputs:** `output/figures/RQ1–4/`, `output/tables/RQ1–4/`, `output/figures/survival/`, `output/tables/survival/`

---

## Finding 1 – CAF2-AR+ depletion is the strongest compositional signal across clinical axes

**Strength: High confidence (FDR-corrected, n=476 ROIs)**

`stromal-CAF2-AR+` is significantly reduced with higher Gleason grade (Spearman rho=-0.29, FDR=4e-7),
in cribriform ROIs (rank-biserial rbc=-0.38, FDR=5e-5), and in stromogenic ROIs.
At the patient level, CAF2-AR+ proportion correlates negatively with max Gleason grade (rho=-0.14, p=0.046).

CAF2-AR+ cells represent smooth-muscle-like or AR-expressing fibroblasts. Their depletion in high-grade
and cribriform disease suggests progressive loss of a "normal" stromal differentiation program.
This is consistent with the draft paper's identification of CAF2 as contractile/SMC-like states,
and with the Bonollo 2020 review describing CAF-marker program acquisition at the expense of
smooth-muscle features in reactive stroma.

Cribriform architecture is a grade-adverse feature with high prognostic relevance. The finding that
CAF2-AR+ loss is the strongest compositional signal for cribriform ROIs (more so than for Gleason grade
alone) points to an interaction between cribriform epithelial architecture and stromal state.

**Uncertainty:** Whether CAF2-AR+ loss is causal or merely tracks with tumor grade requires
experimental validation. Multiple ROIs per patient contribute independently.

---

## Finding 2 – CAF1-CD105+ cells are systematically periglandular; CAF2-AR+ cells are not

**Strength: High confidence (n=475 ROIs with both subtypes, Wilcoxon p=1.9e-75)**

Measured as mean nearest-neighbor distance to luminal epithelium (in 32-pixel units):
- CAF1-CD105+: mean 20.6 px, median 15.0 px  
- CAF2-AR+: mean 41.3 px, median 28.8 px  

Permutation-based contact enrichment (observed/expected at radius 32):
- CAF1-CD105+: 1.17 (slightly above chance)
- CAF2-AR+: 0.49 (strongly below chance; contact with luminal epithelium is LESS than expected by chance)

Periglandular fraction (fraction of CAF cells within 32px of luminal epithelium):
- CAF1-CD105+: 86.4%
- CAF2-AR+: 63.4%

This provides direct quantitative evidence for the "periglandular CAF1 niche" described in the draft
paper. CAF2-AR+ cells are systematically excluded from glandular proximity — consistent with a
model where they occupy interstitial/periductal positions rather than periglandular niches.

The contact enrichment of CAF1+ with luminal epithelium (1.17) being modestly above chance is
important: while CAF1+ cells are spatially close to glands, the contact is not dramatically over-represented
relative to their overall density. This suggests proximity is partly determined by tissue organization
(glands take up a large fraction of the ROI) rather than active recruitment.

**Uncertainty:** Radius-32 neighborhood (~32 μm) may not optimally capture periglandular vs
interstitial distinction. Distance analysis assumes uniform cell density, which may not hold.

---

## Finding 3 – CAF1-CD105+ periglandular fraction predicts clinical progression at the patient level

**Strength: Moderate (p=0.012, n=190 patients, 31 events — survives permutation but borderline given multiple testing)**

Univariate Cox PH: CAF1-CD105+ periglandular fraction vs clinical progression  
HR = 0.717 (95% CI: 0.554–0.929), p = 0.012

**Direction is counterintuitive relative to the draft paper's claim** that the periglandular CAF1 niche
has the "strongest adverse clinical association." Here, *higher* periglandular CAF1 fraction → *lower*
hazard of clinical progression.

Possible explanations:
1. **Metric difference**: The draft paper identifies a discrete spatial niche (cluster of CAF1 cells
   tightly organized around glands). Our metric measures the fraction of all CAF1 cells within 32px
   of luminal epithelium — a continuous measure. High scores could reflect maintained glandular
   architecture even in CAF1-enriched stroma.
2. **Confounding by tissue organization**: In higher-grade tumors (worse outcome), glandular
   architecture is disrupted — fewer recognizable glands, so fewer CAF1 cells can be "periglandular"
   even if the total CAF1 count is unchanged or increased.
3. **CAF1 composition**: The CAF1 category includes both CD105+ and CD105- variants. The periglandular
   CD105+ subpopulation identified by the draft paper as adverse may be specifically the densely-packed
   cluster configuration, not just proximity per se.

KM stratified by median CAF1+ periglandular fraction: p=0.07 for PSA recurrence (trend, high fraction = better).

**Uncertainty:** Low event count (31 for clinical progression), no correction for multiple testing across
scores, confounding by Gleason and other variables not modeled.

---

## Finding 4 – B cells form immune clusters in inflamed ROIs (TLS-like mode), but scatter through tissue in non-inflamed ROIs

**Strength: Moderate-high (neighborhood composition analysis, n=351 ROIs with ≥1 B cell)**

B-cell neighborhood composition differs markedly by inflammation status:
- In **inflamed ROIs**, the local neighborhood of B cells contains:
  - Other B cells: 15.5% (vs 4.7% in non-inflamed)
  - T helper cells: 12.8% (vs 7.9%)
  - T-helper-B-cell mixed states: 8.3% (vs 2.1%)
- In **non-inflamed ROIs**, B cells are primarily surrounded by:
  - Luminal epithelium: 15.7% (vs 7.3% in inflamed)
  - CAF1-CD105-: 10.1% (vs 6.0%)

This suggests two distinct organizational modes:
1. **Inflamed**: B cells cluster with other B cells and T helper cells → TLS-like immune aggregates
2. **Non-inflamed**: B cells scatter through the tissue among epithelial and stromal cells → diffuse infiltration

The TLS-like score (fraction of B cells with ≥1 T cell neighbor) is high across all ROIs (mean 0.814)
because T cells are abundant and B cells tend to co-localize with them in tissue. The difference in
*B-B* colocalization between inflamed and non-inflamed ROIs is more informative.

B-cell proportion (from RQ1: rbc=0.60, FDR=1.6e-15) and immune-mixed (rbc=0.57, FDR=5.8e-14) are
the strongest markers of histological inflammation.

TLS score does not associate with PSA recurrence (p=0.72) or clinical progression (p=0.34) at the
patient level. B-cell proportion also does not survive multiple-testing correction vs survival endpoints.

**Uncertainty:** TLS identification without TLS-specific markers (CXCL13, CCL19) relies solely on
spatial co-localization, which is a limited proxy. The absence of a TLS-outcome association may reflect
limited power (31 clinical events) or genuinely null biology.

---

## Finding 5 – Tumor epithelial composition is the most stable within-patient cell type (ICC=0.64); lymphatic endothelial is the most variable (ICC=0.09)

**Strength: High confidence (ICC analysis, n=154 patients with ≥2 ROIs)**

ICC ranked (descending = more stable):
1. `epithelial-tumor`: ICC=0.64 — the tumor epithelial phenotype is highly reproducible across cores from the same patient
2. `epithelial-luminal`: ICC=0.48
3. `stromal-CAF1-CD105-`: ICC=0.47
4. `undefined`: ICC=0.44
5. `immune-myeloid`: ICC=0.41
...
17. `endothelial-lymphatic`: ICC=0.09 — highly variable within patients

Implication: A single biopsy is informative for the tumor epithelial state (high ICC) but may be
insufficient for immune subtypes like B cells (ICC=0.28) or lymphatic vasculature (ICC=0.09).

Intra-patient JSD (composition divergence between ROIs) does not correlate with Gleason grade
(rho=0.06, p=0.46). Heterogeneity within patients is not systematically higher in high-grade disease.

**Uncertainty:** ICC interpretation depends on group size (k_mean here is ~3 ROIs per patient).
Results may shift with more ROIs per patient.

---

## Summary Table

| Finding | Effect | Evidence | Confidence |
|---------|--------|----------|-----------|
| CAF2-AR+ depleted in high-Gleason | rho=-0.29, FDR=4e-7 (Gleason); rbc=-0.38 (cribriform) | RQ1, n=476 ROIs | High |
| CAF1-CD105+ periglandular vs CAF2-AR+ non-periglandular | 20.6 vs 41.3 px NN distance, p=1.9e-75 | RQ4, n=475 ROIs | High |
| CAF1+ periglandular fraction → better clinical progression | HR=0.72, p=0.012 | Survival, n=190 | Moderate |
| B cells cluster (TLS-like) in inflamed, scatter in non-inflamed | B-B neighborhood 15.5% vs 4.7% | RQ3, n=351 | Moderate-High |
| Tumor epithelial phenotype stable across ROIs (ICC=0.64) | ICC=0.64 vs ICC=0.09 for lymphatics | RQ2, n=154 | High |
| Inflammation richest association category | 13/17 meta-labels significant, FDR<0.1 | RQ1, n=476 | High |

---

## Limitations

- **Multiple ROIs per patient** not accounted for in ROI-level statistical tests (inflates n)
- **Clinical event rates** are modest (54 PSA events, 31 clinical progression events in 190 patients)
  → limited power for survival analyses; findings should be interpreted as hypothesis-generating
- **TLS-like score** relies on spatial proximity, not confirmed TLS markers
- **CAF1+ periglandular finding** has a counterintuitive direction vs the draft paper's niche-based claim;
  requires reconciliation with niche-level analysis and additional confounding control
- All analyses are on tumor ROIs only; non-tumor ROIs were excluded

---

## Recommended Next Steps

1. **Niche-level analysis**: Reconstruct the 18 niches using k-means on local composition vectors
   (as in the draft paper) and directly replicate the association with outcome for the CAF1+ periglandular niche.
   This would resolve the apparent discrepancy with Finding 3.

2. **Confounder adjustment**: Run multivariable Cox models including Gleason grade, age, PSA as
   covariates alongside spatial scores.

3. **Shannon diversity vs clinical variables**: ATHENA-based cell-type Shannon diversity per ROI
   correlated with Gleason and outcome (directly using the starter code in `main.py`).

4. **B-B colocalization density as TLS metric**: Replace TLS score with per-ROI B-B contact density
   or B-cell cluster size, which is more discriminant between inflamed and non-inflamed ROIs than
   fraction of B cells with T cell neighbor.

5. **Image-level representation learning**: Apply CANVAS-like approach (segmentation-free tiles)
   to the IMC images to discover spatial patterns not captured by cell-type composition alone.

---

## Finding 6 – Higher spatial Shannon diversity predicts worse PSA recurrence

**Strength: Moderate (Cox HR=1.38, p=0.040, n=190 patients; 21 ROI-level FDR<0.1 associations)**

Global Shannon diversity (ROI-level) and mean local Shannon diversity (cell-level, radius-32 neighborhood)
are both significantly elevated in inflamed ROIs (rbc≈0.42, FDR=1.6e-7) and in PSA-progressing
patients (local Shannon rbc=0.17, FDR=0.014).

Patient-level Cox PH: global Shannon HR=1.38 (1.01–1.87), p=0.040 for PSA recurrence.
Direction: **higher diversity → higher hazard.** This suggests that spatially mixed
microenvironments — where immune, stromal, and epithelial cells are interleaved — are
associated with aggressive disease. The signal may be partly driven by inflammatory
infiltration (inflamed ROIs are more diverse and in patients with worse outcomes).

Cribriform ROIs show elevated SD of local diversity (rbc=0.40, FDR=8.5e-6) but
reduced epithelial neighborhood diversity (rbc=-0.29, FDR=2.6e-3) — cribriform
epithelium forms tight homogeneous clusters despite heterogeneous surrounding tissue.

Global Shannon is not well-correlated with global composition (it captures spatial mixing
orthogonal to proportions), making it a useful complementary metric.

**Uncertainty:** Shannon diversity and inflammation are correlated; survival signal may be
mediated by inflammation rather than spatial organization per se.

---

## Finding 7 – Paired within-patient analysis: CAF2 (both AR+ and AR-) depleted in high-grade cores; lymphatic endothelium also depleted

**Strength: High confidence (paired Wilcoxon, n=89 matched patients)**

Comparing high- vs low-Gleason TMA cores from the same patient (controlling for patient confounders):
- `stromal-CAF2-AR+`: depleted in high-grade cores (rbc=0.29, FDR=0.075)
- `stromal-CAF2-AR-`: depleted in high-grade cores (rbc=0.34, FDR=0.063)
- `endothelial-lymphatic`: depleted in high-grade cores (rbc=0.33, FDR=0.063)
- `epithelial-tumor`: enriched in high-grade cores (rbc=0.31, FDR=0.070)
- `stromal-other`: enriched in high-grade cores (rbc=0.27, FDR=0.089)

The depletion of BOTH CAF2-AR+ and CAF2-AR- in high-grade cores confirms that the CAF2
depletion is not specific to the AR expression state but reflects a general loss of the
smooth-muscle/contractile CAF2 program in high-grade disease. Lymphatic endothelial depletion
in high-grade cores suggests reduced lymphangiogenesis or loss of organized lymphatic
networks at sites of dedifferentiation.

Spatial Shannon diversity does NOT differ significantly between matched high/low cores
within patients (h_global p=0.15), despite differing between patients with different outcomes.

---

## Finding 8 – The PCa TME spatial wiring diagram: pericyte-vessel contact strongest enrichment; B cell-luminal epithelium strongest exclusion

**Strength: High confidence (476 ROIs, 30 permutations per ROI)**

Mean pairwise contact enrichment matrix (log2 observed/expected) reveals:

*Strongest spatial associations (off-diagonal enrichments):*
- `stromal-pericytes ↔ endothelial-blood-vessel`: +1.24 — pericytes wrap vessels (biological validation)
- `immune-T-cells ↔ immune-mixed`: +1.19 — T cells nest in mixed immune clusters
- `immune-myeloid ↔ immune-mixed`: +0.97

*Strongest spatial exclusions (off-diagonal depletions):*
- `immune-B-cells ↔ epithelial-luminal`: -1.29 — B cells systematically excluded from luminal epithelium
- `epithelial-tumor ↔ stromal-CAF2-AR+`: -1.25 — CAF2-AR+ avoids tumor epithelium (complements RQ4)
- `epithelial-luminal ↔ endothelial-lymphatic`: -1.23

*Self-contact enrichments (spatial clustering tendency):*
- Pericytes (2.18) and immune-mixed (2.03) form the tightest clusters
- Epithelial-luminal (0.77) has the lowest self-contact — most interleaved with neighbors

*Gleason-correlated interaction shifts (30 significant, FDR < 0.05):*
- `CAF1-CD105+ ↔ epithelial-luminal` contact DECREASES with Gleason (rho=-0.23, FDR=1.4e-5)
  → Loss of periglandular CAF1 organization with grade progression
- `CAF2-AR+ ↔ CAF2-AR+` self-contact INCREASES with Gleason (rho=+0.20, FDR=1.9e-4)
  → Remaining CAF2-AR+ cells cluster more tightly as their abundance falls
- `CAF1-CD105+ ↔ epithelial-tumor` also DECREASES with Gleason (rho=-0.19, FDR=3e-4)

*Stromogenic interactions (97 significant):*
- `epithelial-luminal self-contact` strongly DECREASES in stromogenic ROIs (rbc=-0.49)
  → Reactive stroma breaks up luminal epithelial clustering
- `CAF2-AR+ ↔ epithelial-luminal` INCREASES in stromogenic ROIs (+0.44)
  → In reactive stroma, CAF2-AR+ cells are brought closer to epithelium

---

## Finding 9 – Most TME variance is between patients (>70% for major cell types); clustering is continuous not discrete

**Strength: High confidence (variance decomposition, 190 patients)**

Fraction of total cell type variance explained by between-patient differences:
- CAF1-CD105+: 76%, epithelial-tumor: 74%, B cells: 72%, epithelial-luminal: 71%
- Lowest: immune-mixed (56%), immune-fibrocytes (56%), pericytes (59%)

Patient-level hierarchical clustering: silhouette scores are low (max 0.12 at k=2),
indicating **gradual/continuous inter-patient variation** rather than sharp TME subtypes.
Nevertheless all 17 cell types differ significantly across k=3 clusters (KW FDR < 0.1),
confirming real biological structure in the patient-level composition space.

The dominant axis of inter-patient variation (PC1) likely reflects the
immune-hot vs immune-cold and stromogenic vs non-stromogenic dimensions.

---

## Finding 10 – K-means niche reconstruction (k=18): CAF2-depleted niches track Gleason; no periglandular CAF1 niche with adverse outcome

**Strength: High confidence for Gleason/inflammation associations (n=476 ROIs); Null for survival (no FDR-significant Cox results)**

K-means clustering (k=18) on per-cell radius-32 local composition vectors recovered 18 biologically interpretable spatial niches:

*Key niches identified:*
- **N7** (81% epithelial-tumor): pure tumor niche — increases with Gleason (rho=−0.10, FDR=0.082)
- **N17** (74% CAF2-AR+), **N3** (45% CAF2-AR+): CAF2-dominant stromal niches — depleted in high-grade ROIs (N3 rho=−0.25, FDR=8.6e-7; N17 rho=−0.22, FDR=6.7e-6)
- **N13** (36% luminal + 27% CAF2-AR+): mixed luminal-CAF2 niche — depleted in high-grade (rho=−0.19, FDR=1.8e-4) and cribriform ROIs (rbc=−0.29, FDR=8.9e-3)
- **N11** (39% T-cells, 12% immune-mixed): immune niche — strongly enriched in inflamed ROIs (rbc=−0.68, FDR=1.2e-19)
- **N2** (90% luminal): pure luminal gland niche — enriched in non-inflamed, non-cribriform ROIs
- **N6** (35% luminal + 24% CAF1-CD105+): putative periglandular CAF1 niche

*Reconciling with RQ4:*
The periglandular CAF1 niche (N6) does not survive as an adverse survival predictor in Cox PH
(psa_recurrence HR=not computed; no niche reached FDR<0.1). The niche-level analysis is
consistent with RQ4's direction: neither the continuous periglandular fraction nor the discrete
niche usage predicts worse outcome. The discrepancy with the draft paper's claim likely reflects
that the draft identified a *specifically dense/organized* periglandular architecture (using
their 18-niche pipeline tuned for cluster compactness), not just compositional co-occurrence.

*17 significant ROI-level niche associations (FDR<0.1):*
- CAF2-rich niches (N3, N17) decrease with Gleason and in cribriform — replicates Finding 1
- Immune niche (N11) is the single strongest Gleason-adjacent signal for inflammation
- Pure luminal niche (N2, N12) increases in non-inflamed ROIs

*Survival:* No niche Cox results survived FDR correction for either PSA recurrence or clinical
progression. This is consistent with limited power (54/31 events, 190 patients).

**Uncertainty:** MiniBatchKMeans on a 300k-cell subsample may not perfectly replicate the
draft paper's clustering (which used full data + specific seeding). The k=18 choice matches
the paper but was not independently validated here via silhouette/gap statistic.

---

## Finding 11 – B-cell spatial cluster size is a far better TLS proxy than fraction-with-T-cell; nominally associates with PSA recurrence

**Strength: High for inflammation discrimination; Low for survival (p=0.076, not FDR-significant)**

B-cell connected-component clustering in radius-32 graphs (RQ10) resolves RQ3's TLS metric failure:

*ROI-level inflammation discrimination (FDR<0.1):*
- Max B-cell cluster size: rbc=−0.63, FDR=6.4e-18 (inflamed vs non-inflamed)
- B-B contact density: rbc=−0.63, FDR=9.9e-18
- Fraction of B cells in large clusters (≥5): rbc=−0.49, FDR=1.5e-18
- Mean cluster size: rbc=−0.59, FDR=6.9e-16

These metrics are far more discriminative for inflamed ROIs than the RQ3 TLS score (fraction
of B cells with ≥1 T cell neighbor, which was ~0.81 across all ROIs).

*Correlation with RQ3 TLS score:* B-B contact density correlates with the RQ3 TLS score
(Spearman rho=0.447, p=8.7e-25), confirming the old metric captures some of the same signal
but is poorly calibrated.

*Gleason:* None of the B-cell cluster metrics associate significantly with Gleason grade after
FDR correction. B-cell clustering is a signal of inflammation, not grade per se.

*Patient-level survival:* Max B-cell cluster size per patient → PSA recurrence HR=1.002
(95% CI: 1.000–1.004), p=0.076 (direction: larger clusters → higher hazard). This does not
survive FDR correction but is directionally consistent with inflamed TME predicting worse
outcomes (see Finding 6: Shannon diversity also predicted worse PSA recurrence).

**Uncertainty:** The TLS concept requires CXCL13/CCL19 expression and specific architecture
that spatial clustering alone cannot confirm. The cluster size–survival association is weak
and may be mediated by inflammation.

---

## Finding 12 – Shannon spatial diversity and CAF1 periglandular fraction predict PSA recurrence and clinical progression independently of Gleason grade

**Strength: Moderate — Gleason-adjusted Cox, n=190 patients; Shannon p=0.036; CAF1 periglandular p=0.008**

Multivariable Cox models (TME score + Gleason grade; RQ11) reveal which spatial signals
survive adjustment for the strongest clinical confounder:

*PSA recurrence (54 events):*
| Score | HR per SD | 95% CI | p (Gleason-adjusted) |
|-------|-----------|--------|----------------------|
| Global Shannon diversity | 1.38 | 1.02–1.85 | **0.036** |
| Local Shannon diversity | 1.31 | 0.99–1.73 | 0.057 |
| CAF1 periglandular fraction | 0.81 | 0.64–1.02 | 0.077 |
| B-cell max cluster size | 1.20 | 0.98–1.47 | 0.079 |
| CAF2-AR+ proportion | 1.03 | 0.79–1.35 | 0.805 |

*Clinical progression (31 events):*
| Score | HR per SD | 95% CI | p (Gleason-adjusted) |
|-------|-----------|--------|----------------------|
| CAF1 periglandular fraction | 0.69 | 0.52–0.90 | **0.008** |
| Global Shannon diversity | 1.27 | 0.88–1.84 | 0.209 |

**Key interpretation:**

1. **Shannon diversity (global ROI-level) is the strongest Gleason-independent predictor of PSA
   recurrence** (HR=1.38, p=0.036). Higher spatial mixing → worse outcome, even after accounting
   for grade. This suggests spatial disorganization of the TME carries prognostic information
   beyond what is captured by Gleason scoring alone.

2. **CAF1 periglandular fraction is a significant Gleason-independent predictor of clinical
   progression** (HR=0.69, p=0.008). Higher fraction → lower hazard of clinical progression,
   even after Gleason adjustment. This confirms Finding 3 is not simply a Gleason surrogate.
   The mechanistic interpretation remains: high-grade disease loses glandular architecture,
   reducing measurable CAF1 periglandular fraction — the metric acts as a proxy for preserved
   tissue organization. Whether this is causal or organizational requires further investigation.

3. **CAF2-AR+ proportion does not predict survival independent of Gleason** (p=0.80). Its
   cross-sectional association with Gleason (Finding 1) is fully explained by grade.

4. The full model with all 5 scores + Gleason was underpowered (54 events / 6 covariates for
   PSA recurrence, 31 / 6 for clinical progression) and was not fitted.

**Uncertainty:** No FDR correction for multiple tests across scores. Sample size limits power
for independent validation. These are exploratory results requiring prospective validation.

---

## Finding 13 – PC1 of the patient TME space = immune/B-cell axis; PC2 = CAF1 spatial axis; PCs do not predict survival

**Strength: High for biological interpretation; Null for survival**

PCA on patient-level CLR-transformed composition (17 cell types, 190 patients):
- PC1 (24.3% variance): strongest loadings are B cells, immune-mixed, immune-T cells (+) vs
  luminal epithelium, CAF1-CD105+ (−). Interpretation: immune infiltration axis.
- PC2 (17.0% variance): dominated by spatial and stromal features.
- Top 5 PCs explain 68% of total inter-patient variance.

PC1 correlates strongly with B-cell cluster size (rho=−0.57, p=4.7e-18) and Shannon diversity
(rho=−0.45, p=8.6e-11) — confirming it captures the immune-hot/cold dimension. It does NOT
correlate with Gleason grade (rho=−0.11, p=0.13), confirming immune infiltration is largely
independent of grade at the patient level.

No PC predicts PSA recurrence or clinical progression in univariate Cox (all p>0.15).
The dominant axis of inter-patient variation (immune hot/cold) is not prognostic on its own.

---

## Finding 14 – High-grade disease shows CONVERGING inter-patient TME variance (rho=−0.90, p=0.037)

**Strength: Moderate (only 5 grade groups, small n at GG5=10; Levene FDR not significant)**

Between-patient SD of CLR-transformed cell type proportions decreases monotonically with
Gleason grade group (Spearman rho=−0.90 between mean SD and grade, p=0.037 across 5 grade groups).

This supports the **convergence hypothesis**: as tumors progress to higher grade, the TME
converges toward a similar composition — high tumor epithelium, depleted CAF2, reduced
immune infiltration. Low-grade disease (GG1) is more TME-heterogeneous than high-grade (GG5).

Per-cell-type: the `endothelial-lymphatic` variance decreases most strongly with grade
(rho=−1.0 in SD vs grade Spearman, though only 5 points). `undefined` and `stromal-other`
also converge. No cell type shows diverging variance.

**Uncertainty:** Only 10 patients in GG5; Levene's test has limited power per grade group.
No single cell type reaches FDR<0.1. The convergence signal is in the mean across all cell
types (rho=−0.90) rather than any specific cell type.

---

## Finding 15 – Spatial score redundancy: B-cell metrics, niche pairs, and Shannon metrics are highly co-linear

**Strength: High (Spearman correlations, n=190 patients)**

Co-variation matrix of 20 patient-level spatial scores reveals tight redundancy groups:

*Redundant pairs (rho > 0.85):*
- bcell_max_cluster ↔ bcell_bb_density: rho=0.945 — same biology, one is sufficient
- niche_7 ↔ niche_14: rho=0.935 — both tumor-epithelial niches (81% and 49% tumor-epi)
- niche_2 ↔ niche_12: rho=0.876 — both pure luminal niches
- h_global ↔ h_local_mean: rho=0.863 — global and local Shannon capture same variation
- caf1_periglandular ↔ caf1_nn_dist: rho=−0.938 — inverse measures of same proximity

*Strongest cross-group negatives (biological signal):*
- niche_11 (immune) ↔ PC1: rho=−0.70 — confirms PC1 is the immune axis
- h_global ↔ niche_2 (pure luminal): rho=−0.67 — high spatial diversity ↔ few cells
  in homogeneous luminal niches; diversity and architecture are anti-correlated
- bcell_bb_density ↔ PC1: rho=−0.60 — B-cell clustering drives PC1

After FDR correction across 60 score × clinical-variable pairs, 0 pairs survive.
Individual scores that were significant in prior analyses (h_global, caf1_periglandular)
remain nominally associated but multiple testing is too strict with 60 simultaneous tests.

Practical implication: future models should use one representative per redundancy group
(e.g., h_global, caf1_periglandular, caf2_ar_plus_prop, bcell_max_cluster as the minimal set).

---

## Finding 16 – Interaction-based patient clustering is weaker than composition-based (silhouette 0.05 vs 0.12)

**Strength: High confidence (hierarchical clustering, 190 patients, 289 interaction features)**

Patient-level hierarchical clustering on the 17×17 interaction enrichment profiles (RQ14):
- Best silhouette score: 0.050 at k=2 — *lower* than composition-based clustering (0.12, RQ7)
- Spatial wiring profiles are even more continuously distributed than cell-type compositions
- k=2 clusters differ modestly in inflammation status (inflamed: 25 vs 11 patients per cluster)
  but not in Gleason, stromogenic status, or survival (Cox p=0.17 for PSA recurrence)

**Interpretation:** The spatial interaction matrix encodes how cells are wired together, but
this wiring is largely determined by composition (which cell types are present). The interaction
profiles do not contain sufficient patient-discriminating information beyond what composition
already captures. The continuous variation structure seen in RQ7 persists — and is even more
extreme — in the high-dimensional interaction space.

The absence of discrete clusters in both composition and interaction spaces consistently points
to a patient TME landscape that is a continuum, not a set of discrete subtypes.

---

## Finding 17 – CD8+ T cells are rarely spatially excluded; FoxP3+ Tregs track with inflammation and Gleason

**Strength: Moderate (n=476 ROIs; immune exclusion phenotype too rare for powered analysis)**

Immune phenotype classification (based on CD8+ proportion and fraction in epithelial contact):
- Inflamed: 343/476 ROIs (72%)
- Desert: 119/476 ROIs (25%)
- Excluded: 14/476 ROIs (3%)

The "excluded" phenotype — CD8+ cells present but physically barred from tumour epithelium — is very
uncommon in this cohort. This differs from bladder/breast cancer where exclusion is a dominant phenotype.

FoxP3+ regulatory T cells (Tregs) are the strongest spatial predictor among immune subtypes:
- FoxP3+ proportion strongly correlates with inflammation status (rbc=0.62, FDR<1e-10)
- FoxP3+ proportion increases with Gleason grade (rho=0.22, FDR=0.002)
- Univariate Cox: FoxP3 HR numerically extreme and unstable (HR=6e+28, p=0.99) — model is not
  reliably estimable; likely separation or near-zero event overlap with Treg-high ROIs

CD8+ spatial metrics (fraction in epithelial contact, exclusion score) do NOT predict survival
after Gleason adjustment (all p>0.1). Immune composition is less prognostic than stromal architecture
in this dataset.

**Clinical implication:** Immune exclusion is too rare to serve as a stratification biomarker. Treg
enrichment may have utility as a Gleason-independent inflammation marker, but Cox instability
requires a larger cohort for reliable estimation.

---

## Finding 18 – Shannon diversity signal is ~half inflammatory, ~half architectural

**Strength: High (Spearman rho, n=476 ROIs)**

Decomposing the Finding 6 Shannon diversity signal (HR=1.38 for PSA recurrence) by restricting entropy
to non-immune cells:

- Global vs non-immune Shannon: rho=0.939 — they are nearly collinear
- Global Shannon vs immune fraction: rho=0.681 — 46% of Shannon variance explained by immune infiltrate
- Non-immune Shannon Cox (Gleason-adjusted, PSA recurrence): HR≈1.2, p=0.082 — attenuated from p=0.036
- Stromal-only Shannon Cox: p=0.352 — no signal

**Interpretation:** The Shannon diversity prognostic signal is mixed:
~46% is captured by immune infiltration (which drives diversity by adding heterogeneous immune cells),
and ~54% reflects genuine architectural disorganisation of the non-immune tissue compartment.
Removing immune cells weakens but does not eliminate the survival association, suggesting a real
(though modest) architectural component. Shannon is not a pure proxy for either process alone.

This explains why Shannon is partially redundant with immune metrics (B-cell clustering, PC1) in
the score co-variation analysis (Finding 15).

---

## Finding 19 – Dense periglandular CAF1 niches are structurally distinct from diffuse CAF1 proximity

**Strength: Moderate (DBSCAN-based, n=476 ROIs; survival association does not survive Gleason adjustment)**

DBSCAN (eps=40, min_samples=4) on periglandular CAF1-CD105+ cells identifies dense niche clusters:
- Mean fraction of CAF1 cells in dense niche: 26.6%
- `periglandular_frac` vs `frac_caf1_in_dense_niche`: rho=0.075, p=0.10 — essentially uncorrelated

Key clinical associations:
- `periglandular_frac` ↑ with inflammation (rbc=0.22, FDR=0.029) — continuous proximity is
  partly inflammatory co-localisation
- `max_dense_cluster_size` ↑ with Gleason grade (rho=0.12, FDR=0.041) — denser CAF1 clusters
  appear in higher-grade disease (potentially reactive)
- `periglandular_frac` ↓ in cribriform ROIs (rbc=-0.17, FDR=0.088) — confirms Finding 2 direction

Survival (Gleason-adjusted Cox):
- `max_dense_cluster_size`: HR=1.003 per unit, p=0.070 for clinical progression (nominally adverse)
- `periglandular_frac`: HR=0.20, p=0.12 for clinical progression (direction: protective, consistent
  with Finding 3 but not significant after adjustment)
- All dense niche scores: p>0.28 after Gleason adjustment

**Interpretation:** Dense periglandular CAF1 clusters and diffuse periglandular proximity are
structurally distinct (low correlation, opposite Gleason associations), confirming they measure
different biology. The continuous periglandular fraction (Finding 3: HR=0.69) may reflect normal
gland-associated stroma, while dense DBSCAN clusters (slightly adverse Gleason association) may
reflect reactive or desmoplastic responses. Neither metric achieves Gleason-independent survival
significance in this cohort size.

---

## Finding 20 – A spatial signature distinguishes cribriform from non-cribriform ROIs with AUC=0.73

**Strength: High for classification; moderate for individual features (n=476 ROIs, 56 cribriform)**

Univariate associations with cribriform label (Mann-Whitney, FDR-corrected):

| Feature | rbc | FDR |
|---------|-----|-----|
| CAF1–luminal contact enrichment | +0.44 | 7.7e-7 |
| CAF2-AR+ proportion | +0.34 | 1.2e-4 |
| CAF2–luminal contact enrichment | +0.26 | 4.4e-3 |
| Luminal self-contact enrichment | −0.16 | 0.12 (ns) |

Note: positive rbc means the feature is higher in non-cribriform ROIs (i.e., these features are
depleted in cribriform).

Multi-feature logistic regression (5-fold cross-validated AUC):
- Spatial + composition: **AUC = 0.725 ± 0.066**
- Composition only (CAF proportions + epithelial proportion): **AUC = 0.668 ± 0.040**
- **Delta AUC = +0.057** from adding spatial features

Top classifier features by coefficient magnitude:
1. CAF2-AR+ proportion (depleted in cribriform) — strongest single feature
2. CAF1–luminal contact enrichment (lower in cribriform — less organised periglandular niche)
3. Luminal NN distance mean (higher in cribriform — more dispersed gland cells)
4. CAF2–luminal contact enrichment (lower in cribriform)

**Clinical implication:** A computational spatial signature combining stromal depletion markers and
spatial organisation metrics can detect cribriform architecture with AUC=0.73. Spatial features add
meaningful discriminative power (+5.7% AUC) beyond composition alone. This supports the potential
for automated cribriform detection from IMC or spatially-resolved transcriptomics data, complementing
subjective pathologist assessment. The CAF1 spatial organisation of the periglandular niche is a
new structural correlate of cribriform architecture not previously described.
