# Cancer-associated fibroblast classification in single-cell and spatial proteomics data

- Source PDF: `resources/Cords et al. - 2023 - Cancer-associated fibroblast classification in single-cell and spatial proteomics data.pdf`
- Citation: Cords et al., 2023, *Nature Communications* 14:4294
- DOI: `10.1038/s41467-023-39762-1`

## Data and metadata

- Disease / indication: Primary breast cancer discovery cohort with validation across NSCLC, colon cancer, PDAC, and HNSCC datasets.
- Tissue / specimen type: Human breast tumors for discovery and IMC validation; external public stromal / fibroblast datasets from multiple cancers.
- Spatial or single-cell technology: scRNA-seq for discovery; imaging mass cytometry (IMC) for spatial validation.
- Marker panel / assay details: 41-plex IMC antibody panel at 1 micron resolution; stroma-rich and TLS-containing regions were selected.
- Number of patients / samples / regions / images: 14 breast cancer patients in the scRNA-seq discovery cohort; 12 matched breast tumors by IMC; 7-13 IMC regions per patient.
- Number of cells / objects analyzed: About 119,000 total discovery cells including 16,704 stromal cells, of which 14,315 were CAFs and 2,389 pericytes; IMC quantified 553,121 objects.
- Cohorts / validation sets: External public cohorts in lung, colon, pancreas, and head-and-neck cancer; validation included both transcriptomic and protein-level analyses.

## Biological key findings

- The paper defines a CAF taxonomy with 9 CAF phenotypes plus pericytes.
- CAFs separate into broad FAP-positive activated / myofibroblast-like states and FAP-negative structural / perivascular states.
- mCAFs were linked to matrix remodeling, TGF-beta signaling, KRAS signaling, and EMT programs.
- iCAFs carried inflammatory and complement-related programs, including CXCL12, CXCL14, IL6, and IL6-JAK-STAT3 signaling.
- tCAFs and hsp_tCAFs showed stress, hypoxia, and tumor-like programs.
- vCAFs localized near endothelial cells and vessel-like structures, while rCAFs localized near TLS-like immune aggregates.
- ifnCAFs and CD10-positive / CD73-positive tCAFs were especially close to tumor cells.
- Similar CAF phenotypes could be recovered across multiple cancer types, suggesting a reusable cross-cancer classification.

## Methods used to obtain the findings

- Experimental methods: Previously generated breast cancer scRNA-seq, new IMC on matched tumors, and validation using public scRNA-seq cohorts.
- Preprocessing / segmentation / annotation: CellRanger and Seurat for RNA; sctransform, PCA, UMAP, and hierarchical clustering; IMC segmentation with ilastik plus CellProfiler and arcsinh transformation.
- Statistical / computational analysis: Differential expression with MAST, GSEA with `singleseq`, Seurat anchor integration, Harmony, Rphenoannoy, FlowSOM in CATALYST, and differential abundance with edgeR / diffcyt.
- Spatial analysis: Neighborhood analysis within 30 microns, permutation-based interaction scoring, and distance-to-tumor-border analysis.
- Validation: External multi-cancer dataset validation and protein-level confirmation with IMC.

## Code and data availability

- GitHub / code:
  - `https://github.com/BodenmillerGroup/CAFclassification`
  - `https://doi.org/10.5281/zenodo.7540622`
- Data:
  - `https://doi.org/10.5281/zenodo.5769017`
  - `https://zenodo.org/record/7540604`
  - `https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10607/`
  - `GSE132465`
  - `GSE154778`
  - `GSE212966`
  - `GSE103322`
  - `E-MTAB-6149`
  - `E-MTAB-6653`

## Uncertainties / missing metadata

- The total number of IMC regions across all patients was not given explicitly.
- Validation cohorts are heterogeneous and were not summarized with one uniform patient count in the paper text reviewed.
- The IMC panel could not cleanly separate vCAFs from pericytes because RGS5 staining was unavailable.
- hsp_tCAFs were identified in scRNA-seq but not resolvable by the IMC marker panel.

```yaml
machine_readable:
  title: "Cancer-associated fibroblast classification in single-cell and spatial proteomics data"
  year: 2023
  disease: "breast cancer; cross-cancer validation in NSCLC, colon cancer, PDAC, and HNSCC"
  technology:
    - "scRNA-seq"
    - "imaging mass cytometry"
  n_patients: 14
  n_samples: 14
  n_cells: "Discovery: 16704 stromal cells (14315 CAFs); IMC: 553121 objects"
  code_links:
    - "https://github.com/BodenmillerGroup/CAFclassification"
    - "https://doi.org/10.5281/zenodo.7540622"
  data_links:
    - "https://doi.org/10.5281/zenodo.5769017"
    - "https://zenodo.org/record/7540604"
    - "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10607/"
    - "GSE132465"
    - "GSE154778"
    - "GSE212966"
    - "GSE103322"
    - "E-MTAB-6149"
    - "E-MTAB-6653"
```
