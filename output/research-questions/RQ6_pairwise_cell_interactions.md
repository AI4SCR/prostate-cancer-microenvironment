# RQ6: Full Pairwise Cell-Cell Interaction Profiles

## Status: COMPLETED

## Biological Question

Which cell type pairs are spatially enriched or depleted in contact across all 34 cell types?
Do interaction profiles differ between Gleason grades, stromogenic and inflamed states,
and between PSA-recurrent and non-recurrent patients?

## Rationale

RQ4 focused on CAF–epithelial interactions. A full interaction matrix (34×34 cell types,
observed vs permutation-expected contacts) gives a complete spatial wiring diagram of the TME.

Key questions:
- Which interactions are most consistent across patients (low inter-patient variance)?
- Which interactions differ most between high- and low-Gleason ROIs?
- Is there a "spatial interaction signature" that predicts outcome?

## Methods

For each ROI, build radius-32 graph. For each pair of cell types (A, B):
- Count observed A→B contacts
- Permute labels N=50 times, compute expected contacts
- Enrichment score = log2(observed/expected + pseudocount)

Aggregate to patient level (mean or max across ROIs).

## Analysis Plan

### Step 1 — Per-ROI interaction enrichment matrix
- 17×17 meta-label interaction matrix per ROI (tractable; 34×34 is large)
- Permutation-based enrichment (50 permutations per ROI)

### Step 2 — Mean interaction matrix across all ROIs
- Visualize as heatmap — the "average PCa TME wiring diagram"

### Step 3 — Differential interaction analysis
- For each interaction pair: Mann-Whitney between Gleason GG1-2 vs GG3-5
- Between stromogenic vs non-stromogenic
- FDR correction

### Step 4 — Interaction-based patient clustering
- Patient-level interaction vectors → hierarchical clustering
- Are clusters associated with outcome?

## Completed TODOs
- [ ] Per-ROI interaction matrix computed
- [ ] Mean interaction heatmap
- [ ] Differential interactions by clinical group
- [ ] Patient clustering
