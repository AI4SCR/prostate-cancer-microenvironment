# Immune organization in sentinel lymph nodes of melanoma patients is prognostic of distant metastases

- Source PDF: `resources/Amitay et al. - 2024 - Immune organization in sentinel lymph nodes of melanoma patients is prognostic of distant metastases.pdf`
- Citation: Amitay et al., 2024, bioRxiv preprint
- DOI: `10.1101/2024.11.24.625041`

## Data and metadata

- Disease / indication: Melanoma, with a focus on prognosis from sentinel lymph node organization.
- Tissue / specimen type: Retrospective FFPE sentinel lymph node biopsies from a pathology archive; 2-3 annotated ROIs per case were used for TMAs.
- Spatial or single-cell technology: MIBI-TOF spatial proteomics and NanoString CosMx spatial transcriptomics.
- Marker panel / assay details: 39-antibody MIBI panel and 960-gene CosMx panel; MIBI FOVs were 400 x 400 or 800 x 800 microns; CosMx acquisition used 5 fluorescent antibodies.
- Number of patients / samples / regions / images: 69 patients total; 39 negative lymph nodes and 30 metastatic lymph nodes; 68 patients profiled by MIBI, 33 by CosMx; 177 cohort MIBI images plus 16 controls; 77 CosMx images.
- Number of cells / objects analyzed: 1,652,274 MIBI-segmented cells and 848,727 CosMx-segmented cells.
- Cohorts / validation sets: Retrospective four-group cohort stratified by nodal status and distant metastasis outcome; antibody validation on control tissues; CellSighter training / validation split; leave-one-out predictive-model evaluation.

## Biological key findings

- Metastatic sentinel lymph nodes were enriched for myeloid cells, neutrophils, and exhausted cytotoxic T-cell states.
- Patients with nodal metastases who did not later develop distant metastases had more activated and stem-like CD8 T-cell states plus TCF-positive and CD103-positive dendritic-cell programs in the metastatic region.
- Patients who later developed distant metastases showed more CD163-positive / CD209-positive M2-like macrophages and lower tumor-cell HLA-I and SOX10.
- Outside the metastasis, good-outcome metastatic cases had more Tregs in follicles, more effector-memory CD8 T cells in the T-zone, and more PD-L1-positive / CD69-positive myeloid cells.
- In negative sentinel lymph nodes, poor-outcome patients showed more CCR7-positive CD4 T cells, CCR7-positive dendritic cells, Tregs, and PD-L1-positive dendritic cells in the T-zone.
- Good-outcome negative nodes showed sinus expansion, CD206-positive sinus macrophages, plasmablast-sinus proximity, and mast-cell enrichment, consistent with an extrafollicular-like response.
- Low-plex feature models predicted future distant metastases with high AUC in both metastatic and negative-node settings.

## Methods used to obtain the findings

- Experimental methods: Retrospective FFPE cohort, TMA construction, 39-antibody MIBI-TOF profiling, and 960-gene CosMx profiling.
- Preprocessing / segmentation / annotation: MAUI processing with manual denoising; Mesmer segmentation for MIBI; CellPose segmentation for CosMx; CellSighter and InSituType cell typing; manual curation.
- Statistical / computational analysis: T-cell subclustering with k-means, differential expression in Scanpy, Mann-Whitney tests, PCA, linear SVM models, leave-one-out validation, and ROC/AUC evaluation.
- Spatial analysis: DINO embeddings, PhenoGraph clustering of local neighborhoods, logistic regression niche assignment, and macroenvironment mapping across follicles, T-zones, sinuses, and borders.
- Validation: Antibody validation, training/validation image splits, cross-modal comparison between protein and RNA data, and held-out patient prediction.

## Code and data availability

- GitHub / code: Not reported in the paper text reviewed.
- Data: Not reported in the paper text reviewed.

## Uncertainties / missing metadata

- The PDF text reviewed did not include a formal code or data availability section.
- The paper reports modality-specific cell counts, not one deduplicated total.
- Some group-level sample details were visible only in figures or supplements and were not exhaustively extracted here.
- This is a preprint, so journal metadata may change.

```yaml
machine_readable:
  title: "Immune organization in sentinel lymph nodes of melanoma patients is prognostic of distant metastases"
  year: 2024
  disease: "melanoma"
  technology: "MIBI-TOF spatial proteomics and CosMx spatial transcriptomics"
  n_patients: 69
  n_samples: 69
  n_cells: 2501001
  code_links: []
  data_links: []
```
