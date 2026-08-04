# Open questions

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
