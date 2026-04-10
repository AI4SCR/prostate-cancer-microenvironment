# A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging

- Source PDF: `resources/A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging.pdf`
- Citation: Keren et al., 2018, *Cell* 174:1373-1387
- DOI: `10.1016/j.cell.2018.08.039`

## Data and metadata

- Disease / indication: Triple-negative breast cancer (TNBC).
- Tissue / specimen type: Archival FFPE TNBC tissue microarray cores from biopsies; validation and control tissues included tonsil, lymph node, placenta, gastrointestinal tract, skin, breast, and breast carcinoma.
- Spatial or single-cell technology: Multiplexed ion beam imaging by time-of-flight (MIBI-TOF).
- Marker panel / assay details: 36-protein antibody panel; cohort imaging at 2048 x 2048 pixels and 500 nm resolution; platform described as supporting automated 40-plex imaging.
- Number of patients / samples / regions / images: 41 TNBC patients; one patient contributed two cores; 1,763 images total. Patients 1-14 also had an additional H&E core for intra-patient variability assessment.
- Number of cells / objects analyzed: About 3,000-10,000 cells per patient were reported; the paper did not state one cohort-wide total.
- Cohorts / validation sets: Multi-tissue assay validation, serial lymph node sections with split panels, H&E TIL scoring for the full cohort, and duplicate H&E cores for patients 1-14.

## Biological key findings

- Immune infiltration varied by more than two orders of magnitude across TNBC patients.
- Immune composition was structured rather than random, with coordinated co-occurrence of B cells, CD4 T cells, CD8 T cells, NK cells, macrophages, and Tregs.
- Tumors fell into cold, mixed, and compartmentalized spatial archetypes, even when overall immune-cell abundance was similar.
- Spatial enrichment analysis identified conserved patterns including tertiary lymphoid structures, same-lineage clustering, IDO-positive immune clusters, and B-cell depletion at tumor borders.
- PD-1, LAG3, PD-L1, and IDO showed patient-specific and cell-type-specific expression programs.
- Mixed tumors tended to show PD-1 on CD8 T cells with PD-L1/IDO on tumor cells, whereas compartmentalized tumors more often showed PD-1 on CD4 T cells with PD-L1/IDO on immune cells.
- Ordered immune-regulatory structure at the tumor border was associated with improved overall survival.

## Methods used to obtain the findings

- Experimental methods: FFPE TNBC TMA cores were stained with a metal-tagged antibody mix and imaged by custom MIBI-TOF; findings were compared with IHC and H&E pathology review.
- Preprocessing / segmentation / annotation: Spectral calibration, background and artifact removal, DeepCell-based nuclear segmentation, single-cell feature extraction, and FlowSOM-style clustering.
- Statistical / computational analysis: Quantile normalization, k-nearest-neighbor denoising, hierarchical clustering, PCA, chi-square and Wilcoxon tests, randomized null models, Bonferroni correction, and Cox regression.
- Spatial analysis: Pairwise spatial enrichment z-scores, context-aware spatial enrichment, tumor-immune mixing scores, and automated tumor-border detection.
- Validation: Cross-tissue staining specificity checks, serial-section reproducibility, segmentation benchmarking, threshold robustness checks, and H&E-based TIL validation.

## Code and data availability

- GitHub / code:
  - `https://github.com/lkeren/MIBIAnalysis`
  - `http://hub.docker.com/r/vanvalen/deepcell-mibi`
- Data:
  - `https://mibi-share.ionpath.com`

## Uncertainties / missing metadata

- The exact cohort-wide total cell count was not stated.
- Sample count is slightly ambiguous because one patient contributed two cores.
- External links were reported in the paper, but not re-verified here.
- Additional clinical metadata exist in supplementary tables that were not fully extracted.

```yaml
machine_readable:
  title: "A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging"
  year: 2018
  disease: "triple-negative breast cancer"
  technology: "MIBI-TOF spatial proteomics"
  n_patients: 41
  n_samples: 42
  n_cells: "3000-10000 per patient; exact cohort total not stated"
  code_links:
    - "https://github.com/lkeren/MIBIAnalysis"
    - "http://hub.docker.com/r/vanvalen/deepcell-mibi"
  data_links:
    - "https://mibi-share.ionpath.com"
```
