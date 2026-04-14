# RQ2: Intra-patient Heterogeneity Between High- and Low-Grade ROIs

## Status: COMPLETED

## Biological Question

For patients with multiple ROIs sampled from different Gleason areas, how much do ROIs from
the same patient differ in cell type composition and spatial diversity? Is intra-patient
heterogeneity itself prognostically relevant?

## Rationale

The draft paper samples up to four cores per patient from high- and low-Gleason areas. Clinical
sampling of heterogeneous tumors is a major practical challenge. If the TME is highly variable
within a patient, a single biopsy may be misleading. Conversely, if high- and low-grade cores
are concordant in their microenvironment, this suggests TME programming is a patient-level
property.

Key questions:
- Do high-Gleason ROIs from the same patient differ from low-Gleason ROIs in their TME?
- Is intra-patient standard deviation in composition (or Shannon diversity) associated with outcome?
- Are some cell types more variable within patients vs between patients (ICC analysis)?

## Data

- Cell type proportions per ROI (from RQ1)
- Clinical metadata: pat_id, gleason_grp, gleason_pattern_tma_core
- Only patients with ≥2 ROIs at different gleason_grp levels

## Analysis Plan

### Step 1 – Identify multi-ROI patients
- Filter patients with ≥2 ROIs, ideally spanning gleason_grp 1–2 vs 3–5

### Step 2 – Paired comparison per patient
- For each patient: compare cell type proportions between highest and lowest Gleason ROIs
- Wilcoxon signed-rank test across patients

### Step 3 – ICC analysis
- Intraclass correlation coefficient (ICC) per cell type: how much variation is within-patient vs between-patient?
- Cell types with low ICC are highly variable within patients (noisy or locally determined)
- Cell types with high ICC are stable within patients (patient-level program)

### Step 4 – Intra-patient heterogeneity score
- Compute within-patient SD (or pairwise Jensen-Shannon divergence) of ROI composition
- Test whether heterogeneity score associates with outcome

## Expected Outputs

- `output/figures/RQ2/` – scatter plots, ICC bar plots
- `output/tables/RQ2/` – ICC table, heterogeneity scores per patient

## Completed TODOs

- [ ] Multi-ROI patient identification
- [ ] Paired comparison
- [ ] ICC computation
- [ ] Heterogeneity score vs outcome

## Limitations

- "High" and "low" Gleason cores are annotations at the TMA level; may not span the full range
- Small n for some comparisons (patients with multiple grade groups)
