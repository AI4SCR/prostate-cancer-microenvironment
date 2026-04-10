# Characterization of tumour heterogeneity through segmentation-free representation learning on multiplexed imaging data

- Source PDF: `resources/Tan et al. - 2025 - Characterization of tumour heterogeneity through segmentation-free representation learning on multip.pdf`
- Citation: Tan et al., 2025, *Nature Biomedical Engineering*
- DOI: `10.1038/s41551-025-01348-1`

## Data and metadata

- Disease / indication: Lung adenocarcinoma.
- Tissue / specimen type: Patient biopsy cores / TMA IMC images, matched tumor and adjacent normal lung samples, and mouse lung tissue for follow-up experiments.
- Spatial or single-cell technology: Imaging mass cytometry for discovery, with validation by NanoString GeoMx WTA, human and mouse scRNA-seq, and RT-qPCR.
- Marker panel / assay details: 18-channel IMC panel; representative channels include CD20, CD3, CD14, CD68, MPO, HLA-DR, PanCK, DNA1, CD117, CD11c, CD16, CD31, CD4, CD8a, CD94, FoxP3, CD163, and TTF1.
- Number of patients / samples / regions / images: 416 patients in the discovery IMC cohort; 76,963 local tumor microenvironment (LTME) tiles extracted from those biopsies; 143 treatment-naive stage I patients in the validation cohort; 96 GeoMx-profiled tumor cores with 1-3 ROIs each; human scRNA-seq from 18 tumor and 15 normal samples.
- Number of cells / objects analyzed: The main analysis unit was 76,963 LTME objects, each preserving morphology and containing around 20 cells; the total human segmented-cell count was not explicitly stated.
- Cohorts / validation sets: Independent GeoMx/BayesPrism validation cohort, human scRNA-seq validation, mouse KP lung cancer model, and bone-marrow monocyte stimulation with tumor-conditioned media.

## Biological key findings

- CANVAS learned segmentation-free LTME representations directly from multiplexed IMC pixels and grouped them into 50 signatures.
- The resulting signatures captured biologically coherent niches, including B-cell-rich, T-cell-rich, neutrophil-rich, tumor-rich, and morphology-dominant patterns.
- Several signatures were associated with prognosis, smoking status, or stage.
- CANVAS recovered TLS-like immune architectures through co-localized B-cell and T-cell signatures.
- The monocytic signature C10 was the worst-prognosis signature in the discovery cohort.
- C10 was validated in an orthogonal GeoMx/BayesPrism cohort and remained associated with poor progression-free survival.
- Follow-up human and mouse data suggested that tumor-associated monocytes adopt extracellular-matrix-producing programs.

## Methods used to obtain the findings

- Experimental methods: IMC on lung adenocarcinoma biopsies / TMAs, GeoMx WTA on selected ROIs, human tumor-versus-normal scRNA-seq, mouse KP lung model, and monocyte stimulation followed by qPCR.
- Preprocessing / segmentation / annotation: Biopsies were tiled into 64 x 64 micron LTMEs and 5 x 5 micron patches; a masked-image-modeling Vision Transformer was trained with 75% patch masking; segmentation masks were used only for interpretation.
- Statistical / computational analysis: UMAP, clustering into 50 signatures, marker-enrichment scoring, cosine-similarity matching, Kaplan-Meier and log-rank tests, Mann-Whitney U tests with FDR control, BayesPrism deconvolution, and Scanpy / Harmony / Leiden workflows for scRNA-seq.
- Spatial analysis: LTME-to-tissue mapping, graph construction on neighboring LTMEs, adjacency matrices, spectral clustering, and higher-order interaction / homogeneity analysis.
- Validation: Independent GeoMx validation, human tumor-normal scRNA-seq, mouse monocyte profiling, and in vitro tumor-conditioned-media stimulation.

## Code and data availability

- GitHub / code:
  - `https://github.com/tanjimin/CANVAS`
  - `https://doi.org/10.5281/zenodo.14835248`
- Data:
  - `https://zenodo.org/record/7760826`
  - `https://zenodo.org/records/14111081`

## Uncertainties / missing metadata

- The paper did not state one total human segmented-cell count for the discovery IMC cohort.
- The exact number of IMC images or tissue cores per patient was not stated in the extracted text.
- The full 18-marker panel was not listed in one single compact table in the reviewed text.
- The exact total ROI count for the GeoMx validation cohort was not stated.

```yaml
machine_readable:
  title: "Characterization of tumour heterogeneity through segmentation-free representation learning on multiplexed imaging data"
  year: 2025
  disease: "lung adenocarcinoma"
  technology: "imaging mass cytometry with GeoMx, scRNA-seq, and RT-qPCR validation"
  n_patients: 416
  n_samples: 76963
  n_cells: null
  code_links:
    - "https://github.com/tanjimin/CANVAS"
    - "https://doi.org/10.5281/zenodo.14835248"
  data_links:
    - "https://zenodo.org/record/7760826"
    - "https://zenodo.org/records/14111081"
```
