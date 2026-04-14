# RQ11: Multivariable Cox Models — TME Scores Adjusted for Gleason Grade

## Status: COMPLETED

## Biological Question

Do the spatial TME scores identified in earlier analyses (Shannon diversity, CAF periglandular
fraction, CAF2-AR+ proportion, B-cell cluster size) predict clinical outcomes *independently
of Gleason grade*?

## Rationale

All univariate Cox results (RQs 4, 5, 8, 10) may be confounded by Gleason grade, which is
correlated with both TME composition and clinical outcome. Multivariable Cox models with
Gleason as a covariate test whether any TME signal provides independent prognostic information.

## Analysis Plan

### Step 1 — Candidate score list
Collect patient-level scores from prior analyses:
- Shannon diversity (global ROI) — RQ5
- CAF1-CD105+ periglandular fraction — RQ4
- CAF2-AR+ proportion — RQ1
- B-cell max cluster size — RQ10
- CAF1-CD105+ ↔ epithelial-luminal interaction enrichment — RQ6

### Step 2 — Multivariable Cox (Gleason + each score)
- For each score: fit `T ~ E ~ score + gleason_grp`
- Report HR, 95% CI, p-value for score (Gleason-adjusted)

### Step 3 — Full model
- Fit model with Gleason + all scores simultaneously (if n_events allows)

### Step 4 — Report adjusted associations

## Completed TODOs
- [x] Collect patient-level scores (Shannon diversity, CAF1 periglandular, B-cell clusters, CAF2-AR+ proportion)
- [x] Multivariable Cox per score (score + Gleason)
- [x] Full model (skipped — underpowered with 54/31 events for 6 covariates)
- [x] Forest plots and report
