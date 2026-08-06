# Discrepancies: reproduced figures vs. published paper

Visual comparison of `output/figures/*` against the published panels in
`paper.pdf`, done by the user directly. "OK" = visually matches. Items
marked **OPEN QUESTION** need investigation before any fix is attempted --
not yet resolved, not yet guessed at.

## Figure 2

- **2a** (`figure2/figure2a_cell_type_heatmap.pdf`): missing the metagroup
  column annotation present in the paper's version.
- **2b**: OK.

## Figure 3

- **3b**: pericytes are not plotted in yellow (color mismatch vs. paper).
  Not row-clustered, unlike the paper. Overall a different-looking heatmap
  from the published one -- bigger discrepancy than just color/clustering.

## Figure 4

- **4a, 4b, 4c, 4d, 4e**: all OK.

## Figure 5

- **5a**: OK, but the column annotation is placed at the top of the
  heatmap; the paper has it at the bottom.
- **5b**: OK, but the paper's heatmap has no numbers printed in the cells
  (ours does), and the row order differs from the paper.

## Figure 6

- **6a**: **OPEN ISSUE again as of 2026-08-06 -- currently BROKEN.** A fix
  was found and verified on 2026-08-05 (see below), but the script was
  then reverted to the pure verbatim legacy version (per explicit user
  instruction, to diagnose the underlying crash from first principles) and
  was never re-fixed afterward. As currently committed,
  `scripts/figures/figure6_niche_abundance_heatmap.R` crashes with `Error:
  The color mapping should be a named vector or a function.` and cannot
  produce output. The PDF currently sitting at
  `output/figures/figure6/figure6a_niche_proportion_heatmap_tma.pdf` is
  **stale** (left over from the fixed version, before the revert) -- it
  still visually matches the paper, but the current script cannot
  reproduce it.

  Fix found 2026-08-05 (not currently applied): confirmed via direct
  visual comparison (paper page rendered with `pymupdf`) that the paper
  shows no `tma_id` row-annotation track. Root cause of both the missing
  row and the "clustering doesn't look the same" observations at once: an
  earlier fix (uncommenting `tma_id` back into `select()`, to stop a
  genuine 0-row-matrix crash) correctly restored the right row count (459,
  not the 346 you get if `tma_id` is dropped and `distinct()` collapses
  duplicate-clinical-value TMA cores from the same patient) but also left
  `tma_id` visible as a row annotation, which the paper doesn't show. The
  verified fix: keep `tma_id` for indexing only, excluded from the
  displayed annotation data frame. Re-rendered and compared against the
  paper at the time -- row count and overall dendrogram/cluster block
  structure matched closely (same 5 main row clusters, same general niche
  groupings). This fix needs to be re-applied to the current script.
- **6b** (`figure6/niches/km_survival__disease_progr_tumorERG+p53+_ProlifLuminalwith_table.pdf`):
  shows a confidence interval band; the paper's panel does not show a CI.
- **6d**: OK.

## Figure 7

- **7a**: unclear which script should actually produce this panel --
  current mapping in `figures.md`/`table.md` needs re-checking.
- **7b**: **REPRODUCED, RESOLVED 2026-08-05.** Direct visual comparison
  against the published panel (paper page rendered via `pymupdf`) confirms
  the number-at-risk table matches almost exactly: paper shows low
  75/54/20/0, high 115/58/13/0 at t=0/50/100/150; our output
  (`figure7/niches/km_survival__disease_progr_luminal_CAF1(CD105High)with_table.pdf`)
  shows the identical 75/54/20/0 and 115/58/13/0. P-value mismatch
  resolved: the legacy `niche_km.R` (as literally written, un-trimmed
  source) always displays the RAW, unadjusted log-rank p-value on the KM
  plot itself (`ggsurvplot(..., pval = TRUE, ...)`) -- it separately
  computes a BH-adjusted `qval` across ALL niches in one pass (`cols <-
  colnames(df_props)[2:ncol(df_props)]`), but only ever writes that to a
  CSV, never back into the plot. Verified empirically: re-running the KM
  analysis for the full niche family (not just this panel's trimmed 1-2
  niches) reproduces the paper's displayed p-values almost exactly (see
  S6b/S6c below for the confirming numbers). The paper's authors evidently
  substituted the adjusted `qval` into the published panel by hand -- a
  step that exists nowhere in any script. Not a bug in our port; accepted
  as-is per user decision ("as long as the curves match, this is good
  enough").
- **7c** (myCAF): **REPRODUCED, RESOLVED 2026-08-05.** "Not reproduced at all" was
  stale -- predates the stale-save-gate fix made later in the session that
  unblocked this panel; it now saves
  (`figure7/cell_types/km_survival__disease_progr_stromal-CAF1(CD105+)with_table.pdf`).
  Label correctness confirmed directly from the paper's own text (not
  assumed): page 8 states "CD105high are annotated as myCAFs," matching
  our script's `cols <- "stromal-CAF1(CD105+)"` exactly (CD105+ = CD105
  high in this repo's label naming). P-value mismatch resolved by the same
  mechanism as Fig 7b/S6b/S6c below: our raw progression-free p-value is
  0.0137; re-running the full 35-label family (matching legacy
  `celltype_km.R`'s own un-trimmed `cols <- colnames(df_props)[2:ncol(df_props)]`)
  and BH-adjusting jointly gives qval = 0.2389, matching the paper's
  reported p = 0.24 almost exactly.
- **7d, 7e, 7f**: OK.

## Supplementary Figure 1

- **S1a, S1b**: OK.

## Supplementary Figure 2

- **S2a, S2b, S2c**: OK.

## Supplementary Figure 3

- **S3a, S3b, S3c**: OK.

## Supplementary Figure 4

- **S4b**: **INVESTIGATED 2026-08-05, not a data bug.** Spot-checked
  several niches' segment proportions against the paper (rendered via
  `pymupdf`) and they match closely (e.g. niche 1 "luminal" ~85-90%
  epithelial-luminal in both; niche 4 "tumor(ERG+)" ~75-95%
  epithelial-luminal(ERG+) in both). The real difference is orientation:
  the paper's panel b is a horizontal stacked bar (niches as rows), ours is
  vertical (niches as columns) -- harder to compare precisely by eye across
  the two orientations, which is likely why this looked like a proportions
  mismatch. Checked the legacy source directly
  (`000_paper/11_niches/111_heatmaps/visualize_composition.py`): it also
  calls `df_comp_mean.plot(kind='bar', ...)` -- vertical, same as our port.
  So our script is a faithful, verbatim port; the paper's horizontal
  orientation isn't produced by any script found -- almost certainly a
  manual post-processing/rotation step before publication, same category
  as the p-value substitution found for Fig 6b/7b/7c/S6b/S6c. Not fixed
  (would be a disclosed cosmetic deviation from the verbatim `kind='bar'`
  call); flagged here, pending a decision on whether to add it.

## Supplementary Figure 5

- **S5a**: not produced. **Investigated 2026-08-05, reconfirmed genuine
  gap**: re-searched `000_paper/sync_paper/06-spatial-niches/`,
  `000_paper/11_niches/`, and `archive/scripts/11_niches/` (including
  `111_heatmaps/visualize_composition.py`, not previously checked) for any
  script producing the paper's actual S5a panel ("stacked barplot of cell
  type composition of niche 6 ordered by its proportion in each core") --
  found nothing new. `visualize_composition.py` turned out to compute the
  same corpus-wide mean-per-niche composition already ported as
  `figureS4b_niche_mean_composition.py`/`figure6c_niche_composition_filtered.py`,
  not the per-core, niche-6-ordered arrangement S5a's legend describes. No
  legacy script anywhere plots per-core composition this way -- confirmed
  genuinely unimplemented, not an overlooked port. Not attempted, per this
  project's rule against writing new, non-ported plotting code.
- **S5b**: correction -- earlier note in this file that S5b's curves were
  wrong was itself mistaken; user confirmed S5b is correct as produced.

## Supplementary Figure 6

- **S6a**: OK.
- **S6b**: events match the paper; p-value discrepancy **RESOLVED
  2026-08-05**, same finding as Fig 7b above. Empirically confirmed: our
  script's raw displayed p-value for niche 2 (`luminal_infiltrated`,
  progression-free) is 0.393 (matches what was observed, "0.39"); running
  the full 19-niche family and BH-adjusting jointly gives qval = 0.6812,
  matching the paper's reported 0.68 almost exactly. Paper's lack of a CI
  band remains unexplained but accepted (cosmetic, not a data
  discrepancy).
- **S6c**: same mechanism, same full-family qval (0.6812, BH ties at that
  rank) -- consistent with niche 8 (`tumor_CAF1_lymphocytes`) showing the
  same pattern.

## Supplementary Figure 7

- **S7**: **REPRODUCED 2026-08-05.** Source identified earlier this session
  (`sync_paper/06-spatial-niches/construction/kmeans_clustering.py`'s
  ARI-robustness section, `n_runs=50`, `RandomState(42)`, same file as
  Figure 5a's clustering step but a different section), ported to
  `scripts/figures/figureS7_ari_robustness.py`, and run end-to-end
  successfully (~2h15m, 50 k-means(k=24) fits on the 2,051,915-row
  neighborhood graph). No bugs found in this section. Output:
  `output/figures/figureS7/ari_boxplot.png`.
