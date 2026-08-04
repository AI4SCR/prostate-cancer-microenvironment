# Open questions

## Figure 3b CAF heatmap marker list

No legacy script produces the correct 12-marker panel from code alone:

- `000_paper/04_heatmaps/2-cell-types-heatmap.R` (`heatmap.caf()`, this script's original port source): 9 markers.
- `000_paper/04_heatmaps/2-1-cell-types-heatmap.R` (newer sibling, pulled in later): 13 markers, with `pdpn` commented out.
- Published figure (per direct visual confirmation): 12 markers -- the 13-marker list with `pdpn` uncommented and `c_casp3`/`ki_67` removed.

`scripts/figures/figure3_caf_heatmap.R` now hardcodes this 12-marker list directly, since no single legacy script produces it. Flagging as open because the marker set was reconstructed by inference + visual confirmation, not derived from a script we can point to as the source of truth.
