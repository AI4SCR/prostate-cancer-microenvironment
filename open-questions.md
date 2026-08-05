# Open questions

## Audit: strict verbatim compliance across all sync_paper-sourced scripts

Per explicit user instruction: only path/env-var changes are a permitted
deviation from a script's cited legacy source -- everything else (enabling
a commented-out save, dropping a missing column, fixing an undefined
variable, relabeling values) must be reverted, even when that means the
ported script now errors out or produces no/incomplete output. Checked
every script whose `figures.md` citation points at `000_paper/sync_paper/`
for any deviation beyond path changes; reverted 8 found:

- **`figure4_metagroup_barplot.py` / `figure4c_patient_cluster_km.R`**: the
  disclosed "C1-C6 -> P1-P6" relabeling (display-only in 4b, underlying
  value in 4c) doesn't exist in `stacked-frequencies-label.py`/
  `patient_risk_group_km.R` at all -- removed. Cluster labels are now
  literally "C1".."C6" again, not "P1".."P6".
- **`figure6e_immune_risk_score_km.R`**: was still built against the older
  `000_paper/11_niches/113_survival/inflammation_outcome.R`, not the
  currently-cited `000_paper/sync_paper/03_survival/inflammatory_niche_risk_group.R`
  -- these genuinely differ (different raw input, `cell_annotations`/
  `cell_annotation.parquet` vs `clusters_annotated_v2.parquet`; no trailing
  diagnostic CSV; `ggsave()` already active in the source, not commented
  out). Rewritten as a real verbatim port of the current source. Confirmed
  it now fails at `compute_label_frequency()` -- `metadata.parquet` (this
  repo's closest analog to `cell_annotation.parquet`) has no `niche`
  column, matching the cited source's own architecture.
- **`figureS1b_gleason_concordance.R` / `figureS3b_cluster_concordance.R`**:
  reverted the `save_dir`-undefined and stray `p <-` assignment fixes.
  Both now fail with `object 'save_dir' not found`, matching
  `patient-core-heterogeneity.R` exactly.
- **`figureS3c_progression_km.R`**: removed the enabled `ggsave()` and the
  `p_prog_combined` (plot/table) object -- `patient_risk_group_km.R`'s
  progression-free panel only ever does bare `p_prog$plot` (auto-print,
  discarded) then a commented-out `ggsave()` referencing undefined
  `result_dir`. Now fails with `object 'result_dir' not found`.
- **`figureS5b_inflammation_km.R` / `figureS6a_stromogenic_km.R`**: were
  still built against the older `risk_groups_label.R`, not the
  currently-cited `000_paper/sync_paper/03_survival/clinical_risk_group_km.R`
  -- 3 of 4 KM panels there save a combined `plot/table` object, not the
  plot alone, plus different filenames/subdirectory structure. Rewritten as
  real verbatim ports. Both now fail at `ggsave()` with "Cannot find
  directory" -- the source never creates the `dir_inflam`/`dir_stromo`
  subdirectory it saves into, a genuine bug in the current source.

**Deliberately not touched, different category**: Figure 6b/7b/S6bc/7c's
"still implements the older legacy-dir script's logic despite now citing a
`sync_paper` sibling" gap (see the policy audit entry below) is an openly
disclosed citation-vs-implementation mismatch already tracked as a pending
decision, not a silently-applied fix -- out of scope for this pass.
Package substitutions required by genuinely uninstallable dependencies
(`introdataviz` -> `geom_violin`+`position_dodge` in the two violin
scripts) are also not in scope here; those remain under the project's
original, separate "disclosed package substitution" permitted category.

## Supplementary Fig 4b: stacked bars only approximately match the published figure

Per direct user observation, `figureS4b_niche_mean_composition.py`'s output
(mean cell-type composition stacked bar, one bar per niche, all 18 niches)
only approximately matches the published S4b panel -- not confirmed exactly
right, not confirmed wrong in a specific identifiable way either. The
script itself is a verbatim port of `niche_composition.py`'s "MEAN
COMPOSITION BARPLOT" section (see that script's docstring and
`figures.md`'s S4b row) and passed the same fidelity checks as every other
ported script in this project, so this isn't a known logic deviation --
flagging so it isn't silently assumed correct, and as a starting point if
it needs a closer line-by-line recheck (same treatment as the Figure
5a/5b entries below, which are confirmed script-vs-paper mismatches with a
known cause; this one doesn't have a known cause yet). `figure6c_niche_composition_filtered.py`
shares the exact same computation (filtered to niches 6/16/17/18), so
whatever the cause turns out to be likely affects that panel too.

## Policy: `figures.md`'s "Ported from" always cites `sync_paper` when a usable sibling exists

Prompted by a direct question about whether ported scripts used the newest
available source. Checked every "Ported from" entry in `figures.md` that
pointed outside `000_paper/sync_paper/` for a same-purpose `sync_paper`
sibling, and confirmed via git history (not just file mtime, which only
reflects when the legacy repo was pulled/copied into this environment) that
`sync_paper` genuinely holds the later edit in every case: every `sync_paper`
file checked has a last commit dated 2026-04-23 (one bulk "refactoring
paths" commit across the whole tree), versus 2026-03-04/03-17 for every
`000_paper/11_niches/`(or other legacy-dir) counterpart.

Per direct user instruction, adopted as policy and applied to `figures.md`:
**`Ported from` always cites the `sync_paper` sibling when one exists, full
stop -- the older legacy-dir path is cited only when no `sync_paper`
sibling exists at all.** Reproducing the actual script the paper's
collaborators used takes priority over having a working script: if the
`sync_paper` sibling is broken, it's still the cited source and still the
one ported verbatim (bugs and all), documented and flagged for upstream
fixing rather than quietly swapped for an older script that happens to run
(see Figure 6d/7a below -- the one case where this actually applies). This
changed the `Ported from` citation for: Figure 5a (`z_score_heatmap.R` ->
`z_score_composition_vis.R`; `00_kmeans_clustering.py` ->
`kmeans_clustering.py`), Figure 5b (`niche_pairwise_corrleation.R` ->
`correlation_frequencies.R`), Figure 6b/7b/S6bc
(`niche_kaplan_meier_binary.R` -> `niche_km.R`), Figure 6d/7a
(`inflammation_vis.R`/`stromogenic_vis.R`, both now BLOCKED, see below),
Figure 6e (`inflammation_outcome.R` -> `inflammatory_niche_risk_group.R`),
Figure 7c (`label_kaplan_meier_binary.R` -> `celltype_km.R`), and Figure
S5b/S6a (`risk_groups_label.R` -> `clinical_risk_group_km.R`). Figure 6a
already cited its `sync_paper` source (`heatmap_frequencies.R`).

Citing a `sync_paper` source is purely an attribution change for most of
these rows -- the sibling was confirmed functionally identical (or the
script logic was already unaffected) for Figure 5a, 5b, 6e, and S5b/S6a, see
the checks below. **Figure 6b/7b/S6bc and 7c have since been rewritten as
genuine verbatim ports of the sources they cite** (see below) -- this
surfaced that those sources are themselves broken against this repo's
data, so all four panels are now BLOCKED rather than producing (possibly
subtly wrong) output.

**Figure 7c / 6b / 7b / S6bc -- RESOLVED (now genuinely verbatim, but BLOCKED at runtime).**
All four scripts were rewritten as real verbatim ports of their cited
`sync_paper` sources (`celltype_km.R`, `niche_km.R`), trimmed to the
relevant label/niche(s) only. This surfaced why the older implementation
had drifted in the first place -- the current sources are not simply
"the same logic plus BH q-values," they're broken as literally written
against this repo's data:
- `niche_km.R` reads `paths$cell_annotations` (this repo's `metadata.parquet`),
  which has no `niche` column (that only exists in `clusters_annotated_v2.parquet`);
  fails inside `compute_label_frequency()` before reaching any KM code. All
  three niche scripts (6b/7b/S6bc) fail here identically.
- Separately, every plot-saving call in `niche_km.R` is commented out
  regardless (`#print(p_os)`, `#print(p_prog)`, and even its one
  conditional `ggsave()` for specific niches) -- only two aggregate CSVs
  with BH-adjusted q-values would ever be written, and only if the script
  got past the missing-column failure above.
- `celltype_km.R` has its own, unrelated bug: `dir.create(figures, ...)`
  references an undefined `figures` variable (the defined one is
  `figures_dir`) -- fails immediately, before touching any data.
  `figure7c_myCAF_km.R` now fails here, matching exactly.

So the original threshold/BH-adjustment questions (75th vs 50th percentile,
raw vs FDR-corrected p-values) are moot for now -- none of these four
scripts can produce any output at all as genuinely-verbatim ports of their
current sources. Revisit if the underlying data gap (a `niche`-bearing
per-cell table under `EXPORT_DIR`, or a fixed `dir.create` typo upstream)
is ever resolved.

**Checked and confirmed identical / no logic change needed:**
- Figure 6e: sibling `03_survival/inflammatory_niche_risk_group.R` is
  functionally identical (same threshold, same risk_group logic) and
  independently applies the same two disclosed fixes already made in our
  port -- enabling the commented-out `ggsave()` calls, and dropping the
  trailing diagnostic CSV that selects the nonexistent `clinical$sample_name`
  column (see the Figure 6e entry below). This cross-check increases
  confidence in both fixes.
- Figure S5b / S6a: sibling `03_survival/clinical_risk_group_km.R` covers the
  same two sections as `risk_groups_label.R`'s inflammation/stromogenic parts
  (its patient-cluster-KM section was split out separately into
  `patient_risk_group_km.R`, already our Figure 4c/S3c source) with
  `ggsave()` already enabled and no other logic change.
- Figure 5b: sibling `06-spatial-niches/abundance/correlation_frequencies.R`
  confirmed identical, see that entry below.
- Figure 5a: sibling `06-spatial-niches/composition/z_score_composition_vis.R`
  reads the same precomputed `niche_heatmap_data.parquet` our current script
  already reads (not a from-scratch recompute on either side, contrary to an
  earlier pass over this file); the only real differences are cosmetic
  (colors sourced from `resources/colormaps.yaml`-equivalent YAML instead of
  hardcoded, one dead `cell_fun` block dropped, output written once instead
  of twice under two names) -- confirmed via full diff, no logic deviation.
  Its clustering-step sibling, `06-spatial-niches/construction/kmeans_clustering.py`,
  inlines the same `perform_kmeans_clustering()`/`wrapper_nhood_filtering()`
  calls (same `k=24`, `seed=686`, `KMeans(n_init='auto', init='k-means++')`)
  that `figure5_niche_clustering.py` already imports live from the staged
  `utils.clustering` module -- same underlying function, inlined vs.
  imported, not a different implementation.
- Figure 2 heatmap/UMAPs, Figure 4d-e Cox HR, Figure S1a cohort summary: no
  `sync_paper` sibling exists at all (`sync_paper` has no `04_heatmaps`,
  `02_umaps`, or `01_clinical_metadata` equivalent, and no Cox-hazard-ratio
  script under `03_survival`) -- current sources are the only ones available,
  citation unchanged.

## Figure 6d / 7a: BLOCKED -- verbatim ports of the latest (`sync_paper`) source, which errors as literally written

Was: an earlier pass on this question kept `figure6_inflammation_violin.R`/
`figure7a_stromogenic_violin.R` pointed at the older, working
`000_paper/11_niches/112_violinplots/inflammation_vis.R`/`stromogenic_vis.R`,
reasoning that their `sync_paper` siblings were broken and therefore not
"the relevant sibling." **Overridden per direct, explicit user instruction:
reproducing the figures means using the actual script the paper's
collaborators used, verbatim, even if that script is broken -- a broken
script gets documented and escalated to the collaborators to fix, not
silently swapped for a different, older script that happens to run.**

Both scripts are now 1:1 ports of the `sync_paper` source (paths translated
to this repo's `.env` vars, `utils_R.R`'s `compute_label_frequency()`
inlined; no other change), and both are confirmed **BLOCKED** by actually
running them in this environment:

- **Figure 6d** (`000_paper/sync_paper/06-spatial-niches/histology/inflammation_vis.R`):
  fails with `Error: object 'df_long' not found` at
  `df_long$inflammation <- factor(df_long$inflammation, ...)` -- `df_long`
  is referenced but never defined anywhere in the file. Traced via git log:
  its definition was present as of commit `a7d0b602` (2026-04-12 21:03,
  "descriptive niches added") and deleted 12 minutes later by `c185c9e9`
  ("added niches and histology part") while the line using it was left
  behind -- an accidental regression, never fixed by the later `aa8a8b59`
  ("refactoring paths", 2026-04-23) commit, which only touches
  path/config-loading boilerplate.
- **Figure 7a** (`000_paper/sync_paper/06-spatial-niches/histology/stromogenic_vis.R`):
  fails with `Error in c(no = "#9efa70", yes = "#c55797", ) : argument 3 is
  empty` at the `my_cols <- c(...)` definition -- the trailing comma is
  syntactically legal in R (unlike the "unparseable" characterization in
  the now-superseded entry for the older script, see below) but produces an
  empty third argument, a runtime error. A second, independent bug sits
  further down (`p_full <- ggplot(df_long_full, ...)` references
  `df_long_full`, never defined in this file) but isn't reachable/confirmed
  since execution halts at the `my_cols` error first.

Neither bug is fixed here, and the older, working legacy-dir scripts are
not substituted. These are now blocked exactly like the other
missing-data/permission-denied blockers in this file (Figure 6a, Figure
7d-f) -- the fix has to come from whoever maintains `sync_paper`. If it's
fixed upstream, re-pull and re-run; until then these two panels produce no
output. The stale PDFs previously produced by the old (working, but no
longer cited) script version were removed from `output/figures/figure7/`
since they no longer correspond to any current script's output; Figure 6d
had none staged.

Also note: even a fixed version of these scripts would drop the
split-violin/facet-by-up-down-direction panel present in the older scripts,
replacing it with a differently-specified pair of panels (per-niche
`facet_wrap`, then a non-split dodged violin) -- a real panel-scope
difference, not just a bug, worth knowing about if/when this unblocks.

**Not deeply audited, flagged for awareness only (no `figures.md` citation exists to change here):**
- Figure 5a's clustering-support sibling, `06-spatial-niches/construction/kmeans_clustering.py`,
  also bundles `compute_spatial_nhood.py`-adjacent functions (GPU/Leiden
  graph construction, nhood-filtering) belonging to an earlier upstream step
  (building the cell-cell neighborhood graph itself) that no script in this
  repo reproduces either way -- not relevant to the citation change above,
  noted in case that upstream step is ever tackled.
- Figure 7d-f circos: correcting an error in an earlier pass over this
  question, which assumed this panel was still blocked -- it's actually
  already resolved (see the Figure 7d-f entry below) and
  `figure7def_circos_plots.py` runs end to end today.
  `000_paper/sync_paper/06-spatial-niches/interactions/interaction_compute_circos.py`
  is a combined, single-file (582-line) version of the *compute* stages
  (`compute_interactions.py`/`visualize_interactions_lfc.py`-equivalent) --
  the two upstream stages our current script deliberately skips, since their
  output already exists precomputed (per the resolved entry). It doesn't
  contain a `visualize_circos_plot.py`-equivalent (the actual
  figure-drawing stage `figure7def_circos_plots.py` ports), so it isn't a
  relevant sibling to switch to; citation unchanged.

## ~~Figure 6a: wrong legacy script ported, correct one's data is missing~~ PARTIALLY RESOLVED

Was: `scripts/figures/figure6_niche_abundance_heatmap.R` ported
`000_paper/11_niches/111_heatmaps/patient_heatmap.R`, which is
**patient-level** (`props_niche_pat_id.parquet`, rows = patients) -- wrong
granularity versus the paper's actual "per tumor core" Fig 6a legend, which
is also why the row annotations didn't match what's published.

Correct source confirmed by direct user identification:
`000_paper/sync_paper/06-spatial-niches/abundance/heatmap_frequencies.R`
(TMA/core-level; row annotations `pat_id, os_status, disease_progr,
gleason_grp, inflammation, stromogenic_smc_loss_reactive_stroma_present`
-- exact match). (`000_paper/11_niches/111_heatmaps/proportion_heatmap.R`
was the other TMA-level sibling candidate, ruled out -- it has two extra
annotation columns, `gs_grp`/`d_amico_risk`, beyond what's actually shown.)

`figure6_niche_abundance_heatmap.R` has been rewritten as a strict verbatim
port of `heatmap_frequencies.R`. Confirmed by running it that its only
issue is the missing input data below -- fails at exactly one point, the
`read_parquet(".../niche_frequencies_per_tma_id.parquet")` call.

**Legacy bug, preserved verbatim per explicit user instruction (no
disclosed fix -- reverted from an earlier version of this port that did
fix it)**: `tma_id` is commented out of `df_metadata`'s `select()` but
still referenced two lines later, `matrix[as.character(df_metadata$tma_id), ]`.
Since that column doesn't exist, this evaluates to `matrix[character(0), ]`
-- a 0-row matrix. So even once the missing data below is available, this
script as literally ported will not produce a populated heatmap; matches
legacy's own behavior exactly, not "fixed."

**Still blocked**: `niche_frequencies_per_tma_id.parquet` does not exist
anywhere accessible -- checked
`/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood` (nothing matching
`niche_frequencies*` at all) and
`/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/frequencies/` (has a
`niche_frequencies_per_tma.parquet`, no `_id`, but uses an older,
pre-"revised annotation" niche naming scheme, e.g. `TLS_Bcells_Tcells`,
`canonical_BLepithelium` -- not the same data, not a usable substitute).

Per `CLAUDE.md`'s top-of-file hard constraint: **we do not write scripts to
compute missing data, only port scripts against data that already
exists.** The script is ported and correct; it will run once
`niche_frequencies_per_tma_id.parquet` (or the `stacked_frequencies.py` run
that would produce it, per that script's own `group_var='tma_id'` branch)
turns up staged somewhere.

## Figure 3b CAF heatmap marker list

No legacy script produces the correct 12-marker panel from code alone:

- `000_paper/04_heatmaps/2-cell-types-heatmap.R` (`heatmap.caf()`, this script's original port source): 9 markers.
- `000_paper/04_heatmaps/2-1-cell-types-heatmap.R` (newer sibling, pulled in later): 13 markers, with `pdpn` commented out.
- Published figure (per direct visual confirmation): 12 markers -- the 13-marker list with `pdpn` uncommented and `c_casp3`/`ki_67` removed.

`scripts/figures/figure3_caf_heatmap.R` now hardcodes this 12-marker list directly, since no single legacy script produces it. Flagging as open because the marker set was reconstructed by inference + visual confirmation, not derived from a script we can point to as the source of truth.

## Figure 3b CAF heatmap color scale

The two legacy `heatmap.caf()` versions also disagree on the color scale:

- `000_paper/04_heatmaps/2-cell-types-heatmap.R` (our current port source): `colorRamp2(c(-2, 0, 2), c("lightblue", "white", "lightcoral"))`.
- `000_paper/04_heatmaps/2-1-cell-types-heatmap.R` (newer sibling): `colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))`.

Our port uses the older `lightblue`/`white`/`lightcoral` scale. Unconfirmed against the published figure which is correct -- unlike the marker-list question above, this hasn't been visually checked yet.

## Pericytes in CAF UMAP visualizations

Pericytes (`stromal-pericytes`) are excluded from `figure3_caf_umap.py`'s population (`label` must contain `"CAF"`), matching the only CAF-specific UMAP script found in the legacy repo (`000_paper/02_umaps/0-umaps-cafs.py`). If the published Figure 3a actually includes pericytes, there must be a modified/different version of this script that does -- an exhaustive search of the legacy repo (including the recently-pulled `sync_paper/` content) found no such variant. Location of any such modified script is unknown; flagging so this isn't silently assumed correct if the published figure turns out to include pericytes.

## Figure 5a column annotation position: top in legacy code, bottom in published figure

Both legacy copies of the generating script (`000_paper/11_niches/111_heatmaps/z_score_heatmap.R`
and the older `05_nhoods/PCA_NHOODs_clean/R_visualization/z_score_heatmap.R`) construct
`col_ha = HeatmapAnnotation(celltype = ..., main_group = ..., ...)` and pass it as
`top_annotation = col_ha` -- identical in both copies, so this isn't a
stale-vs-updated-script disagreement like the Figure 3b marker list. Per direct
user observation, the published figure has this same annotation drawn at the
*bottom* of the heatmap instead.

`scripts/figures/figure5_niche_heatmap.R` matches the legacy source exactly
(`top_annotation`), per the project's verbatim-port mandate -- not changed to
`bottom_annotation` to visually match the paper, since no legacy script
produces that layout. Flagging as an unresolved script-vs-published-figure
discrepancy rather than silently correcting it.

## Figure 5b: numbers and column order deviate from published figure

`scripts/figures/figure5_niche_correlation.R` matches the legacy source
exactly: all three copies found in the legacy repo -- the port's source
`000_paper/11_niches/111_heatmaps/niche_pairwise_corrleation.R`, its older
twin at `05_nhoods/PCA_NHOODs_clean/R_visualization/niche_correlation/niche_pairwise_corrleation.R`,
and `000_paper/sync_paper/06-spatial-niches/abundance/correlation_frequencies.R`
-- build the identical `Heatmap(corr_matrix, ..., cluster_rows = TRUE,
cluster_columns = TRUE, cell_fun = ...)` call: per-cell correlation numbers
drawn via `cell_fun`, row/column order left to `ComplexHeatmap`'s default
clustering.

Per direct user observation, the published Figure 5b heatmap has no numbers
in the cells, and its row/column order runs the opposite direction from what
this call produces (on a symmetric correlation matrix, `cluster_rows`/
`cluster_columns = TRUE` cluster rows and columns identically, so "reversed
column order" means the whole heatmap is mirrored along its diagonal, not
just one axis). Not changed -- per the project's verbatim-port mandate, same
as the Figure 5a annotation-position case above, since no legacy script
produces a no-numbers/mirrored-order version of this heatmap.

## Figure 5c / S5a: no legacy plotting code for per-core composition bars

Two panels remain genuinely unimplemented (not just unported) after reading
the actual published figure legends (`2026.04.30.721907v1.full.pdf`):

- **Figure 5c legend** (verbatim): "Representative IMC cores with their
  corresponding H&E images... Adjacent stacked barplots show the mean cell
  type composition of the dominant niches in these cores (containing at
  least 15 cells of the niche)." A **per-representative-core, x=niche** bar
  -- requires picking specific representative cores (non-code-derived, same
  as Fig 1/3c-e) and then computing per-core composition for just that
  core's dominant niches.
- **S5a legend**: "Stacked barplot of cell type composition of niche 6
  ordered by its proportion in each core." Per-core cell-type composition
  *within* niche 6, cores ordered by niche 6's abundance.

`niche_composition.py`'s `calculate_freqs()` computes the per-sample data
that *would* support either panel (`freq_nhood_celltype_sample`), but no
script anywhere in the legacy repo actually plots it this way -- only the
corpus-wide mean-per-niche composition (`df_comp_mean`), which is what's
portable and is now implemented as `figureS4b_niche_mean_composition.py`
(all 18 niches) and `figure6c_niche_composition_filtered.py` (niches
6/16/17/18, Figure 6c's top panel only). Figure 6c's own "two example
cores" bottom panel has the same problem as 5c/S5a and is also not
implemented. Not attempted, per the project's rule against writing new,
non-ported plotting code to fill visualization gaps.

**This superseded an earlier, wrong identification of Figure 5c's source**:
before reading the actual figure legend, `000_paper/sync_paper/06-spatial-niches/abundance/stacked_frequencies.py`
was flagged as the correct 5c source, since it's the only script in the two
niche-relevant directories that produces a *TMA-level niche-abundance*
stacked bar (x=TMA, per-TMA p53-niche-high filtering) -- structurally
plausible, but the actual published legend describes a *per-core, x=niche*
composition bar instead, which `stacked_frequencies.py` does not produce.
`stacked_frequencies.py` also has a genuine legacy bug (its annotation panel
indexes a `disease_progr` column that `tma_cols` never keeps -- a `KeyError`
as literally written), which no longer matters for this decision but is
worth remembering if it turns out relevant elsewhere.

**Where to look**: per direct user instruction, for any niche-related figure
work, restrict the search for candidate legacy scripts to
`000_paper/11_niches/` and `000_paper/sync_paper/06-spatial-niches/` first
(and their `/users/amarti51/projects/PCa/...` mirror). Scripts elsewhere
(e.g. `000_paper/100_other_visualization/plot_stacked_frequencies.py`, which
has hardcoded `/users/mensmeng/workspace/...` sys.path imports and looks like
an earlier/abandoned draft) are lower-confidence and shouldn't be assumed
correct just because they surface in a text search.

## Niche numbers 1-18 are not stored in any data file (inferred mapping)

The paper's "niche 1".."niche 18" numbering is purely editorial -- it isn't
a column in `niche_annotations_v2.csv` (whose `cluster` column is the raw
0-23 k-means id; multiple raw clusters map to one named niche) or anywhere
else in the data. The mapping used throughout this repo's newer scripts
(`figure7b_niche9_km.R`, `figureS6bc_niche_km.R`) is inferred from the fixed
`niche_order` list used consistently across
`figure5_niche_heatmap.R`/`niche_composition.py`/`z_score_heatmap.R`:
position 1 = `luminal`, ..., position 18 = `TLS` (full list in
`figureS4b_niche_mean_composition.py`'s `NICHE_ORDER`). This is
independently confirmed by textual cross-reference for niche 6
(`tumorERG+p53+_ProlifLuminal`, matches Figure 6b's already-validated panel)
and niches 16-18 (`immune_bloodvessels_CAF1(CD105-)`/
`Macrophages_Tcells_CAF1(CD105-)`/`TLS`, matches `inflammation_outcome.R`'s
`cols_inflamed` and Figure 6e's legend exactly), giving high confidence
overall for the list-position mapping. Niches 2 (`luminal_infiltrated`), 8
(`tumor_CAF1_lymphocytes`), and 9 (`luminal_CAF1(CD105High)`) -- used for
Figure 7b and Supplementary Fig 6b-c -- are consistent with their paper-text
descriptions but are **not independently confirmed** the way 6/16-18 are.
Flagging in case the KM results for these three don't match the published
p-values.

## myCAF identification: `stromal-CAF1(CD105+)`

Used for Figure 7c (`figure7c_myCAF_km.R`). The paper text states "CD105high
are annotated as myCAFs"; `resources/colormaps.yaml`'s `label:` key confirms
`stromal-CAF1(CD105+)` is the CD105-high CAF1 variant. This corrects an
earlier, wrong guess of `stromal-CAF2(AR+)` (an androgen-receptor-positive
variant, unrelated to CD105 status) made before checking the paper text
directly.

## Figure 6e: `clinical$sample_name` doesn't exist in this repo's export

`inflammation_outcome.R`'s trailing diagnostic CSV (`df_histo`, not the
actual KM panel) selects `clinical$sample_name` -- a column this repo's
`clinical.parquet` (via `scripts/data/export.py`) doesn't produce (checked:
`sample_id`, `tma_sample_id`, `unique_tma_sample_id_{1-4}`,
`napari_sample_id` exist, `sample_name` doesn't). Dropped from the select in
`figure6e_immune_risk_score_km.R` since it's incidental to the panel -- both
KM plots save successfully before this block runs, confirmed by running the
script both before and after the fix. Flagging as an open schema question
(does legacy's `sample_name` correspond to one of this repo's `*sample_id*`
columns, or was it dropped from the export entirely?) rather than guessing
which one to substitute.

## ~~Figure 7a: legacy syntax bug (trailing comma)~~ SUPERSEDED (script no longer ports this source)

Was: `stromogenic_vis.R`'s (`000_paper/11_niches/112_violinplots/`)
`my_cols <- c(no = "#9efa70", yes = "#c55797",)` has a trailing comma inside
`c(...)` -- fixed in `figure7a_stromogenic_violin.R` by removing it, on the
claim that trailing commas make R argument lists unparseable. Also noted
there: legacy's `p_up` panel (stromogenic-specific up/down niches) is
computed and printed but never saved via `ggsave`, matched verbatim.

Superseded: `figure7a_stromogenic_violin.R` no longer ports this script (see
"Figure 6d / 7a: BLOCKED" above) -- per explicit user instruction, it now
ports the latest available source, `000_paper/sync_paper/06-spatial-niches/histology/stromogenic_vis.R`,
which has the identical trailing-comma bug, left unfixed there. Actually
running the current script shows the "unparseable" claim above was
imprecise: R parses a trailing comma in `c(...)` fine (unlike some other
languages), it just produces an empty argument, which errors at evaluation
time (`argument 3 is empty`), not parse time -- everything before that line
runs first. Kept as a historical note in case the `11_niches` version's fix
is ever relevant again.

## `patient-core-heterogeneity.R`: two undefined-variable bugs

Both fixed as disclosed, minimal, save-enabling changes (same precedent as
`figure6_niche_abundance_heatmap.R`'s already-enabled commented-out
`pdf()`/`dev.off()`):
- Both the Gleason-concordance and cluster-concordance sections reference
  `save_dir` in their `ggsave()` calls, but only `output_dir`/`figures_dir`
  are ever defined -- fixed by using `figures_dir` (`figureS1b_gleason_concordance.R`,
  `figureS3b_cluster_concordance.R`).
- The Gleason-concordance section's final `ggsave(plot_path, p, ...)`
  references `p`, which is never assigned (only auto-printed via
  `(g_strip + g_main) + plot_layout(...)`) -- fixed by assigning that
  expression to `p` first (`figureS1b_gleason_concordance.R`).

## ~~Figure 7d-f (circos interaction plots): blocked on missing data~~ RESOLVED

Was: `compute_interactions.py` (first stage of the pipeline) requires
per-sample anndata pickles at
`/users/mensmeng/workspace/nhoods/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas/`,
inaccessible (`Permission denied`) and not staged under `LEGACY_DATA_DIR`.

Resolved -- the same anndata pickles exist at a different, accessible live
path, `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood/CellCellNeighborhoods/anndatas/`
(2,515 files; loadable with a disclosed `anndata`-version compatibility
monkeypatch, see below). But per direct user instruction, the two upstream
stages that would consume them (`compute_interactions.py`,
`visualize_interactions_lfc.py`) were skipped entirely, since neither
produces a paper panel itself -- they're intermediate data-generation steps.
Their combined output already exists precomputed at the exact path the
actual figure-producing script (`visualize_circos_plot.py`) reads:
`.../5-niches/visualization/interactions/redo/per_niche_lfc_above_median/dataframes_v2/{niche}.parquet`
(confirmed present and schema-correct for niches 2/8/9, staged into
`LEGACY_DATA_DIR` for self-containment -- see `data/assets.md`). Ported as
`figure7def_circos_plots.py`, bundling niches 2/8/9 (d/e/f) into one
script per the same "figure-level bundling where it's the natural unit"
precedent as `figureS6bc_niche_km.R`. Verified running end to end.

**Two things worth remembering if `compute_interactions.py`/
`visualize_interactions_lfc.py` are ever ported** (not currently needed,
since their output is already available precomputed):
- The `anndata` pickle version-compatibility issue: pickled with an older
  `anndata` whose internal `AnnDataFileManager` used the state key `_adata`;
  the installed `anndata` 0.12.10 expects `_adata_ref` and raises
  `KeyError` otherwise. Confirmed these objects aren't actually file-backed
  (`_filename`/`_file` are `None`), so it's a pure key-rename fix, not a
  real backed-file problem.
- A sample-ID namespace mismatch: `compute_interactions.py` indexes
  `clusters_annotated.parquet` by the anndata pickle's short napari-style
  filename (e.g. `240217_005`), but that file is indexed by the long-form
  `sample_id` used everywhere else in this repo. Bridging via
  `clinical.parquet`'s `napari_sample_id` column works for 515/541 (95%)
  of pickles; the rest have no matching clinical row (presumably excluded
  from the final cohort).
- **Important, initially-missed nuance**: `compute_interactions.py` reads
  `clusters_annotated.parquet` (no `_v2`), which disagrees with the
  verified-correct `_v2` file on 73% of rows' `niche` column -- this looks
  like a bug at first glance, but isn't: `visualize_interactions_lfc.py`
  computes everything using the old (non-`_v2`) niche assignment throughout
  (this is genuinely what produced the paper's results), then at the very
  end remaps *only the output filenames* from old to new niche names via a
  lookup on the raw k-means cluster ID (present in both files), before
  saving to `dataframes_v2/`. Naively substituting `_v2` throughout the
  computation -- which is what "fix the stale annotation" would naturally
  suggest -- would actually be the real deviation, since it would change
  which cells get grouped into which niche during the interaction
  computation itself.

## Supplementary Fig. 7 (50-seed ARI robustness sweep): no plotting script found

Data exists (`.../PCa_NHood/robustness/*_kmeans_robustness.pkl`) but no
script that plots it (boxplots of pairwise adjusted Rand index per run) was
found anywhere in the legacy repo, including the two niche-relevant
directories above. Would need to be written from the raw sweep output if
implemented -- flagged as a gap, not attempted.

## ~~Figure 4c KM styling: add_censor_mark()/add_pvalue() are not in the legacy source~~ RESOLVED

Was: `ggsurvfit` substitution required reconstructing `survminer::ggsurvplot()`'s argument-driven censor marks/p-value display as separate calls, with unconfirmed fidelity. Resolved -- `survminer` now installs successfully in this environment (see REPRODUCIBILITY.md's install note; `Deriv` pinned to 4.2.0 from CRAN's archive). All three scripts that had substituted `ggsurvfit`/`survfit2` for `survminer::ggsurvplot()` (`figure4c_patient_cluster_km.R`, `figure6_niche_abundance_heatmap.R`, `figure6_km_niche6.R`) reverted to the literal legacy calls. No package-substitution deviation remains for KM plotting anywhere in this repo.

## ~~Figure 4a: no legacy script found, reconstructed from Methods text~~ RESOLVED

Was: an earlier version of `figure4_patient_clustering.py` was a Methods-text
reconstruction (`sns.clustermap`) flagged as not derived from any actual
legacy script. Resolved -- the real generating script was found
(`000_paper/sync_paper/05-heterogeneity/stacked-frequencies-label.py`,
`group_var == 'pat_id'` branch: a stacked bar plot with a dendrogram drawn on
top, not a clustered heatmap). `figure4_patient_clustering.py` was rewritten
from scratch as a verbatim port and validated against the published figure.

## ~~Figure 4c: which of two same-named `risk_groups_label.R` files?~~ RESOLVED (moot)

Was: `000_paper/100_other_visualization/risk_groups_label.R` and
`000_paper/11_niches/113_survival/risk_groups_label.R` have the same
filename but different content, and it was unclear which (if either) mapped
to Figure 4c. Resolved -- neither is the actual source. The real Figure 4c
generating script turned out to be
`000_paper/sync_paper/03_survival/patient_risk_group_km.R`, found after the
user provided it directly. `figure4c_patient_cluster_km.R` is a verbatim
port of that script and is validated.

## ~~Figure 2/3/6 logic-fidelity audit (11-script sweep)~~ RESOLVED

Was: an 11-subagent audit found real logic deviations in
8 of 11 `scripts/figures/*` files versus their legacy sources -- most
seriously `figure2_umap.py` (6 major: wrong `n_neighbors`, missing
normalization, missing marker exclusion, extra explicit `metric=`, narrowed
marker-plot scope, unverified data-loading equivalence), `figure2_cell_type_heatmap.R`
(8 major: wrong filter target, missing subsampling, flipped matrix
orientation, hardcoded clustering params, missing annotations, non-deterministic
colors, unconditional aggregation, undisclosed new output), and
`figure3_caf_umap.py` (3 major: wrong cell filter, wrong `n_neighbors`,
missing normalization) -- plus smaller issues in `figure5_niche_clustering.py`,
`figure5_niche_annotation.py`, `figure5_niche_heatmap.R`,
`figure6_niche_abundance_heatmap.R`, and `figure6_inflammation_violin.R`
(missing `stat_summary` layer). Resolved -- all fixed and re-verified against
current script content (`n_neighbors=50`, `normalize(..., exclude_zeros=True)`,
correct marker exclusion/CAF filter, correct heatmap filter/orientation/annotations/colors,
`stat_summary` layers present in all three violin plots). `figure5_niche_correlation.R`,
`figure6_km_niche6.R`, and `figure4de_cox_hazard_ratio.R` (then `figure4_survival.R`)
passed with 0 real discrepancies at audit time and remain unchanged.
