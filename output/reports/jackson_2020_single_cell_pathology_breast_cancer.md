# The single-cell pathology landscape of breast cancer

- Source PDF: `resources/s41586-019-1876-x.pdf`
- Citation: Jackson et al., 2020, *Nature*
- DOI: `10.1038/s41586-019-1876-x`

## Data and metadata

- Disease / indication: Breast cancer.
- Tissue / specimen type: FFPE breast tumor tissue arranged as TMAs; the discovery cohort used one 0.8 mm tumor core per patient, and the validation cohort used four 0.6 mm cores from distinct tumor regions.
- Spatial or single-cell technology: Imaging mass cytometry (IMC).
- Marker panel / assay details: 35-biomarker panel covering ER, PR, HER2, Ki-67, epithelial, mesenchymal, immune, endothelial, signaling, oncogene, and epigenetic targets.
- Number of patients / samples / regions / images: The abstract reports 352 patients and 720 images. The methods describe a discovery cohort of 281 patients with 381 images and a validation cohort of 72 patients with 344 additional images; the reporting summary mentions 71 validation patients, so counts are slightly inconsistent across sections.
- Number of cells / objects analyzed: 855,668 cells in the first cohort and 411,410 in the validation cohort, for 1,267,078 combined cells if summed.
- Cohorts / validation sets: Discovery cohort from University Hospital Basel and an independent validation cohort from University Hospital Zurich.

## Biological key findings

- Single-cell IMC identified 18 single-cell pathology (SCP) subgroups that refined standard clinical breast-cancer classes.
- SCPs were prognostically distinct, separating favorable HR-positive groups from poor-outcome CKlowHRlow, hypoxic, basal, proliferative, and p53-positive / EGFR-positive states.
- Tumor architecture was organized into discrete communities, and greater spatial heterogeneity was associated with poorer outcomes.
- Stromal environments carried independent prognostic information, including a poor-outcome fragmented environment with proliferative vimentin-high fibroblasts.
- Microenvironment communities added prognostic information beyond single-cell phenotypes alone.
- The independent validation cohort reproduced the main metaclusters and subgroup structure.
- Around 40% of tumors had identical subgroup classification across regions, while about 60% showed at least one discordant region, highlighting spatial heterogeneity.

## Methods used to obtain the findings

- Experimental methods: IMC on FFPE sections using a 35-antibody panel; Hyperion acquisition; TMA-based sampling in the discovery cohort and multi-region sampling in the validation cohort.
- Preprocessing / segmentation / annotation: TIFF conversion, pixel classification with Ilastik, segmentation with CellProfiler, single-cell feature extraction with histoCAT / MATLAB, and spillover compensation with CATALYST.
- Statistical / computational analysis: PhenoGraph clustering, hierarchical metaclustering, t-SNE, Kaplan-Meier and Cox models, log-rank tests, Shannon entropy, Kullback-Leibler divergence, and Pearson correlation for cross-cohort matching.
- Spatial analysis: Nearest-neighbor graphs, permutation-based neighborhood enrichment, graph / community detection for tumor and microenvironment communities, and stromal-environment grouping.
- Validation: Independent external cohort, matched analytical pipeline, cross-cohort cluster matching, and multi-region tumor sampling.

## Code and data availability

- GitHub / code:
  - `https://github.com/BodenmillerGroup/SCPathology_publication`
- Data:
  - `https://doi.org/10.5281/zenodo.3518284`

## Uncertainties / missing metadata

- Patient and image counts are slightly inconsistent across the abstract, methods, and reporting summary.
- The combined total of 1,267,078 cells is derived by summing cohort-specific counts, not quoted as one number in the paper text.
- Some analyses excluded patients with too few communities, so not every downstream model used the full cohort.
- Supplementary tables contain additional patient-level metadata that were not fully extracted here.

```yaml
machine_readable:
  title: "The single-cell pathology landscape of breast cancer"
  year: 2020
  disease: "breast cancer"
  technology: "imaging mass cytometry (IMC)"
  n_patients: 352
  n_samples: 720
  n_cells: 1267078
  code_links:
    - "https://github.com/BodenmillerGroup/SCPathology_publication"
  data_links:
    - "https://doi.org/10.5281/zenodo.3518284"
```
