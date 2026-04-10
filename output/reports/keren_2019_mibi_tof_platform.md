# MIBI-TOF: A multiplexed imaging platform relates cellular phenotypes and tissue structure

- Source PDF: `resources/MIBI-TOF A multiplexed imaging platform relates cellular phenotypes and tissue structure.pdf`
- Citation: Keren et al., 2019, *Science Advances*
- DOI: `10.1126/sciadv.aax5851`

## Data and metadata

- Disease / indication: Triple-negative breast cancer for the main biological application; the paper also includes platform demonstrations in hippocampus, GI tract, melanoma, colorectal cancer, and tonsil tissue.
- Tissue / specimen type: Archival human FFPE tissue sections, including a whole TNBC section and additional validation/demo tissues.
- Spatial or single-cell technology: MIBI-TOF with dynamic SIMS and orthogonal TOF detection.
- Marker panel / assay details: 36-plex TNBC panel; 34-plex CRC panel; 40-antibody hippocampus panel; 20-plex GI tract demo; platform capability described up to 42 simultaneous metal-labeled antibodies.
- Number of patients / samples / regions / images: Main TNBC biological analysis used one archival TNBC section from one patient, tiled across 520 400-micron fields of view and rescanned in 8 high-resolution regions.
- Number of cells / objects analyzed: 11,120 cells across 8 TNBC regions; a CRC subanalysis compared 109 beta-catenin-high and 115 beta-catenin-low tumor cells.
- Cohorts / validation sets: Prior 41-patient TNBC cohort referenced for comparison; technical validation included NanoSIMS comparison, foil sensitivity testing, GaAs dynamic-range testing, and absolute reporter quantification.

## Biological key findings

- Tumor-cell phenotypes were region-dependent within a single TNBC section.
- Nearby regions were more similar than distant ones, and regional phenotype distributions deviated from Poisson expectation.
- Immune-cell composition was more stable across regions than tumor-cell composition.
- Regional immune co-occurrence patterns matched the group’s earlier 41-patient TNBC cohort.
- PD-L1 was predominantly found on immune cells, while PD-1 was predominantly found on CD4 T cells across the analyzed regions.
- Spatial analysis identified a concentric tumor-border organization with PD-L1-positive / IDO-positive myeloid-like cells adjacent to tumor cells and lymphocyte-rich zones farther out.
- In CRC, nuclear beta-catenin-high tumor cells were more proliferative and more IDO1-positive than beta-catenin-low cells.

## Methods used to obtain the findings

- Experimental methods: FFPE sectioning, deparaffinization, antigen retrieval, metal-conjugated antibody staining, and MIBI-TOF imaging.
- Preprocessing / segmentation / annotation: Denoising, DeepCell segmentation, arcsinh-transformed per-cell expression normalized by cell area, FlowSOM clustering, and Cytobank-based gating for the CRC subanalysis.
- Statistical / computational analysis: t-SNE, clustering, Pearson correlation of regional composition versus distance, Poisson expectation tests, FDR-corrected KS tests, chi-square tests, and Mann-Whitney tests.
- Spatial analysis: Tumor-immune mixing scores, randomized pairwise spatial-enrichment z-scores, and automated tumor-border identification.
- Validation: Instrument sensitivity benchmarking, stability tests, absolute quantification of conjugated reporters, and comparison against prior TNBC biology.

## Code and data availability

- GitHub / code:
  - `https://github.com/lkeren/MibiAnalysis`
  - `www.deepcell.org` for segmentation access, as noted in the paper
- Data:
  - `https://mibi-share.ionpath.com`

## Uncertainties / missing metadata

- This is partly a platform paper, so counts are reported for representative experiments rather than a single study-wide table.
- The total number of cells in the referenced 41-patient TNBC cohort was not reprinted here.
- Demo datasets are technical validations and not directly comparable with the main TNBC analysis.

```yaml
machine_readable:
  title: "MIBI-TOF: A multiplexed imaging platform relates cellular phenotypes and tissue structure"
  year: 2019
  disease: "triple-negative breast cancer"
  technology: "MIBI-TOF"
  n_patients: 1
  n_samples: 8
  n_cells: 11120
  code_links:
    - "https://github.com/lkeren/MibiAnalysis"
    - "www.deepcell.org"
  data_links:
    - "https://mibi-share.ionpath.com"
```
