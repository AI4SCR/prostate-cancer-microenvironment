# Cross-Paper Methods Overview

This document is a working synthesis of the main methodological approaches used across the papers in `output/reports/`, excluding the internal draft `martinelli_2026_pca_draft.md`.

The goal is not to define one canonical pipeline. The goal is to capture:

- what kinds of methods were used
- what parameter values or ranges were reported
- what normalization and spatial-analysis choices recurred
- what aspects are worth checking in a new dataset

## Papers included

- Keren et al. 2018, TNBC MIBI-TOF
- Keren et al. 2019, MIBI-TOF platform paper
- Jackson et al. 2020, breast cancer IMC
- Danenberg et al. 2022, breast cancer TME structures with IMC
- Cords et al. 2023, CAF classification from scRNA-seq and IMC
- Amitay et al. 2024, melanoma sentinel lymph nodes with MIBI-TOF and CosMx
- Tan et al. 2025, segmentation-free IMC representation learning with GeoMx/scRNA-seq validation
- Bonollo et al. 2020, prostate CAF review

## 1. Study design and cohort scale

### Common study designs

- Single-center retrospective TMA cohorts with linked clinical metadata
- Spatial profiling of FFPE sections or TMAs
- Discovery cohort plus one of:
  - external validation cohort
  - orthogonal modality validation
  - internal serial-section / technical validation
- Cross-cancer transfer or taxonomy validation using public scRNA-seq datasets

### Typical cohort sizes reported

- Small / technical demonstration:
  - 1 patient, 8 high-resolution regions in Keren 2019
- Medium translational cohorts:
  - 41 TNBC patients in Keren 2018
  - 69 melanoma patients in Amitay 2024
  - 14 discovery patients plus 12 IMC validation cases in Cords 2023
- Large breast cancer TMA cohorts:
  - 352 patients, 720 images in Jackson 2020
  - 693 analyzed patients, 749 images from 794 tissue spots in Danenberg 2022
- Large representation-learning discovery cohort:
  - 416 IMC-profiled lung adenocarcinoma patients in Tan 2025

### Typical data volumes

- Tens of thousands of cells:
  - 11,120 cells in Keren 2019
- Hundreds of thousands of objects:
  - 553,121 IMC objects in Cords 2023
  - 855,668 cells in the Jackson discovery cohort
- Millions of cells:
  - 1,652,274 MIBI cells in Amitay 2024
  - 848,727 CosMx cells in Amitay 2024
  - 1,267,078 combined cells in Jackson 2020

### Aspects to investigate in a new dataset

- Is the study unit patient, core, ROI, image, tile, or cell?
- How many patients are truly analyzable after QC?
- Are multiple cores / ROIs per patient available?
- Is there an independent validation cohort or only internal validation?
- Are clinical outcomes mature enough for survival analysis?

## 2. Tissue type, sampling, and imaging unit

### Common tissue formats

- FFPE TMAs were the dominant design in most studies.
- Whole-section imaging appeared in some platform or focused biology studies.
- Some studies mixed TMA-based discovery with orthogonal validation on different platforms.

### Common spatial analysis units

- Single cells after segmentation
- Whole ROIs / fields of view
- Local neighborhoods or patches
- Higher-level graph communities / niches / structures

### Reported ROI / image / patch scales

- 400 x 400 micron or 800 x 800 micron MIBI fields in Amitay 2024
- 520 tiled 400-micron fields plus 8 rescanned regions in Keren 2019
- 1 micron IMC resolution in Cords 2023
- 500 nm image resolution in Keren 2018
- 64 x 64 micron local tumor microenvironments with 5 x 5 micron patches in Tan 2025
- 0.6 mm to 1.0 mm TMA spots in Danenberg 2022
- 0.8 mm discovery cores and four 0.6 mm regional cores in Jackson 2020

### Aspects to investigate

- Is the relevant biological unit a gland, invasive front, core, or microenvironment tile?
- What is the expected spatial scale of the signal:
  - immediate cell-cell contact
  - perivascular neighborhood
  - gland-level niche
  - whole-core composition
- Are the imaged ROIs representative or preselected for interesting biology?

## 3. Assays and modality combinations

### Spatial proteomics platforms

- MIBI-TOF:
  - Keren 2018
  - Keren 2019
  - Amitay 2024
- Imaging mass cytometry:
  - Jackson 2020
  - Danenberg 2022
  - Cords 2023
  - Tan 2025

### Multi-omic or orthogonal validation modalities

- scRNA-seq:
  - Cords 2023
  - Tan 2025
- GeoMx WTA:
  - Tan 2025
- CosMx:
  - Amitay 2024
- H&E or pathology scoring:
  - Keren 2018
  - Jackson 2020
  - Danenberg 2022

### Reported panel sizes

- 18-channel IMC in Tan 2025
- 35-biomarker IMC in Jackson 2020
- 36-protein MIBI panel in Keren 2018
- 36-plex TNBC MIBI panel in Keren 2019
- 37-dimensional IMC in Danenberg 2022
- 39-antibody MIBI panel in Amitay 2024
- 41-plex IMC panel in Cords 2023
- Additional examples in Keren 2019:
  - 34-plex CRC
  - 40-antibody hippocampus
  - 20-plex GI tract

### Working range seen in these papers

- Spatial proteomics panels commonly fell in the 18-41 marker range.
- The most common practical range was roughly 34-39 markers.

### Aspects to investigate

- Does the panel cover:
  - epithelial identity
  - stromal / CAF states
  - endothelial / vascular states
  - major immune lineages
  - functional states such as proliferation, apoptosis, checkpoint biology, oncogenic state
- Are markers sufficient to resolve desired subtypes, or will important states collapse?
- Are niche-defining markers missing, for example:
  - pericyte markers
  - TLS maturation markers
  - macrophage polarization markers
  - fibroblast lineage markers

## 4. Preprocessing, denoising, and segmentation

### Common preprocessing steps

- Spectral calibration and background removal in MIBI workflows
- TIFF conversion before segmentation and downstream processing
- Denoising or artifact removal before cell calling
- Spillover compensation in IMC workflows
- Manual QC and artifact filtering

### Tools reported

- DeepCell:
  - Keren 2018
  - Keren 2019
- Ilastik + CellProfiler:
  - Jackson 2020
  - Cords 2023
  - Danenberg 2022
- Mesmer:
  - Amitay 2024
- CellPose:
  - Amitay 2024 for CosMx
- MAUI with manual denoising:
  - Amitay 2024

### Common segmentation choices

- Nuclear segmentation followed by cell expansion or mask-based cytoplasm assignment
- In CRC analysis in Keren 2019, nuclei were expanded by 10 pixels to define cytoplasmic annuli.
- In Amitay 2024 CosMx, cells were expanded by 25 pixels after segmentation.
- Segmentation-free learning was used in Tan 2025 for the main representation pipeline, with segmentation masks retained only for interpretation.

### Common segmentation caveats

- Missing markers limited separation of nearby stromal states:
  - vCAF versus pericyte in Cords 2023 because RGS5 was unavailable
- Some subtypes identified in scRNA-seq could not be resolved in spatial protein panels:
  - hsp_tCAF in Cords 2023

### Aspects to investigate

- Is nucleus-only segmentation adequate, or do you need membrane-aware segmentation?
- Are cell boundaries biologically reliable in the tissue type?
- Which compartments can be separated confidently with the available markers?
- What failure modes are common:
  - merged adjacent cells
  - stromal cell fragmentation
  - under-segmentation in dense lymphoid regions
  - misassignment at tumor-stroma borders

## 5. Feature extraction and normalization

### Common feature representations

- Mean per-cell marker intensity within segmentation masks
- Area-normalized per-cell expression
- Per-ROI or per-neighborhood composition vectors
- Graph-derived connectivity profiles
- Learned image embeddings from raw multiplexed pixels

### Normalization and transformation choices reported

- Arcsinh transformation:
  - Keren 2019
  - Cords 2023 IMC pipeline
- Quantile normalization:
  - Keren 2018
- 99th percentile censoring to remove outliers:
  - Jackson 2020
- 99.9th percentile clipping:
  - Martinelli draft only, excluded from this synthesis
- sctransform for scRNA-seq:
  - Cords 2023
- Seurat anchor-based integration and Harmony for batch correction:
  - Cords 2023
  - Tan 2025 for scRNA-seq preprocessing
- Per-cell normalization by cell area:
  - Keren 2019

### Practical ranges from these papers

- Outlier clipping or censoring was typically done at high percentiles rather than using aggressive hard thresholds.
- Area normalization appeared when cell size was expected to confound intensity.
- IMC / MIBI pipelines often used transformed protein intensities plus composition summaries rather than raw counts directly.

### Aspects to investigate

- Are intensities comparable across slides, batches, and acquisition days?
- Is area normalization needed for large stromal cells versus compact lymphocytes?
- Should normalization happen:
  - per image
  - per slide
  - per batch
  - globally
- Do outliers represent biology or technical artifacts?

## 6. Cell typing and clustering

### Common cell-typing strategies

- Unsupervised clustering of single-cell marker profiles
- Manual merging of clusters into interpretable phenotypes
- Semi-supervised or supervised cell typing in transcriptomic data
- Cross-modal label transfer from scRNA-seq to spatial proteomics

### Reported tools

- FlowSOM / FlowSOM-style clustering:
  - Keren 2018
  - Keren 2019
  - Cords 2023 via CATALYST
- PhenoGraph:
  - Jackson 2020
  - Danenberg 2022
  - Amitay 2024 for local microenvironments
- UMAP:
  - Cords 2023
  - Tan 2025
- t-SNE:
  - Keren 2019
  - Jackson 2020
- CellSighter:
  - Amitay 2024
- InSituType:
  - Amitay 2024
- Seurat / Harmony:
  - Cords 2023
  - Tan 2025

### Typical outputs

- 10 recurrent TME structures in Danenberg 2022
- 18 SCP subgroups in Jackson 2020
- 9 CAF phenotypes plus pericytes in Cords 2023
- 50 LTME signatures in Tan 2025

### Aspects to investigate

- What is the desired resolution:
  - broad compartments
  - functional subtypes
  - clinically actionable rare states
- Which clusters are reproducible across patients?
- Which clusters are likely technical splits rather than biological states?
- Can the phenotype labels be supported directly by markers present in the panel?

## 7. Spatial analysis frameworks

### Common spatial approaches

- Pairwise neighborhood enrichment / avoidance
- Distance-to-border or distance-to-compartment analysis
- Tumor-immune mixing or compartmentalization scores
- Local neighborhood composition clustering
- Community detection on cell graphs
- Higher-order niche / structure definitions
- Segmentation-free embedding of local patches

### Reported parameter values and design choices

- 30 micron neighborhood radius in Cords 2023
- 4 micron adjacency in Jackson 2020
- Randomized null models for spatial enrichment:
  - Keren 2018
  - Keren 2019
  - Jackson 2020
  - Cords 2023
- Tumor-border analysis:
  - Keren 2018
  - Keren 2019
  - Cords 2023
- DINO embeddings plus PhenoGraph in Amitay 2024
- Local tiles of 64 x 64 microns and patching into 5 x 5 microns in Tan 2025
- Spectral clustering of LTME interaction graphs in Tan 2025

### Common derived spatial outputs

- Cold / mixed / compartmentalized tumor architectures
- Tumor and microenvironment communities
- CAF-rich, immune-rich, or epithelial-rich niches
- TLS-like structures
- Perivascular structures
- Tumor-border myeloid or lymphoid rims

### Aspects to investigate

- Is the biology mostly contact-based, radius-based, border-based, or graph-based?
- What is the right neighborhood radius for the tissue architecture?
- Should niches be defined from:
  - single-cell interactions
  - local composition vectors
  - graph communities
  - image embeddings
- Are niche labels interpretable and stable across parameter changes?

## 8. Statistical modeling and outcome association

### Common association tests

- Mann-Whitney or Wilcoxon rank-sum tests
- Chi-square tests
- Pearson or Spearman correlation
- KS tests
- Generalized linear models
- Differential abundance models

### Survival and clinical modeling

- Kaplan-Meier curves:
  - Keren 2018
  - Jackson 2020
  - Danenberg 2022
  - Tan 2025
- Cox proportional hazards models:
  - Keren 2018
  - Jackson 2020
  - Danenberg 2022
- Logistic regression for subtype or niche prediction:
  - Danenberg 2022
  - Amitay 2024 for niche assignment
- Random forest validation:
  - Danenberg 2022
- Linear SVM classification:
  - Amitay 2024
- Differential expression / abundance tools:
  - MAST
  - edgeR
  - diffcyt

### Multiple-testing control

- Benjamini-Hochberg or FDR correction was common.
- Bonferroni correction appeared in Keren 2018.

### Aspects to investigate

- Which level should be used for inference:
  - cell
  - ROI
  - patient
- How should multiple ROIs per patient be aggregated?
- Are survival models adjusted for grade, subtype, stage, or treatment?
- Are rare-state associations robust after multiple-testing correction?

## 9. Validation strategies

### Common validation patterns

- Independent patient cohort validation:
  - Jackson 2020
  - Danenberg 2022
  - Tan 2025
- Cross-cancer validation:
  - Cords 2023
- Orthogonal modality validation:
  - Amitay 2024 with MIBI and CosMx
  - Tan 2025 with IMC, GeoMx, scRNA-seq, and mouse follow-up
- Technical validation:
  - Keren 2019 instrument and sensitivity benchmarking
  - Keren 2018 staining and segmentation reproducibility
- Pathology-based validation:
  - Keren 2018 H&E TIL scoring

### Common robustness checks

- Serial-section reproducibility
- Randomized null models for spatial interactions
- Repeated clustering or repeated assignment
- Out-of-sample prediction across sites or centers
- Concordance with pathology scoring or known tissue structures

### Aspects to investigate

- Can the main findings reproduce:
  - across centers
  - across slides
  - across modalities
  - across tissue regions
- Are the niches biologically plausible on visual inspection?
- Are the most important findings dependent on one marker, one slide, or one patient subgroup?

## 10. Data and code reporting patterns

### Common reporting patterns in published papers

- GitHub repositories for analysis code
- Zenodo for data or code snapshots
- Public repository accessions for transcriptomic datasets
- Web portals for browsable image data

### What to record in a new study

- Raw-image availability
- Segmentation masks
- Per-cell tables
- ROI metadata
- Clinical metadata at the level allowed by governance
- Code for preprocessing, cell typing, spatial analysis, and figure generation
- Exact versions of external models and tools used for segmentation or embedding

## 11. Method families seen across the papers

### A. Marker-driven single-cell spatial proteomics

Used in:

- Keren 2018
- Keren 2019
- Jackson 2020
- Danenberg 2022

Common shape:

- segment cells
- quantify marker intensities
- cluster cells into phenotypes
- analyze pairwise spatial organization
- relate patterns to clinical outcome

### B. Spatial proteomics plus transcriptomic reference mapping

Used in:

- Cords 2023
- Amitay 2024
- Tan 2025

Common shape:

- learn or define cell states in transcriptomic data
- validate or project them into spatial data
- study spatial proximity or niche membership
- test clinical association

### C. Segmentation-free local representation learning

Used in:

- Tan 2025

Common shape:

- tile images into local windows
- learn embeddings directly from pixels
- cluster local microenvironments
- validate the learned states with orthogonal modalities

### D. Review / conceptual synthesis

Used in:

- Bonollo 2020

Common use:

- identify candidate stromal states, signaling axes, and biological hypotheses
- not suitable as a direct protocol source for numeric defaults

## 12. Provisional defaults and ranges to consider for a new dataset

These are not recommendations yet. They are the main values repeatedly encountered in the papers and worth testing.

### Panel design

- Start by checking whether the panel falls within the observed 18-41 marker range.
- A practical “rich spatial proteomics” range in these papers was about 34-39 markers.

### Spatial unit

- ROI sizes in the hundreds of microns were common.
- If using local neighborhoods, 30 microns is one concrete radius reported.
- If using graph adjacency, a very local 4 micron adjacency was reported in one breast IMC study.

### Normalization

- Test arcsinh-transformed intensities for spatial proteomics.
- Test per-cell area normalization where cell size differs strongly across compartments.
- Evaluate high-percentile censoring or clipping of extreme intensities.

### Clustering

- FlowSOM, PhenoGraph, and graph-based or embedding-based clustering were the dominant approaches.
- UMAP or t-SNE were common for inspection, not as final evidence by themselves.

### Spatial inference

- Use permutation-based nulls for neighborhood enrichment whenever possible.
- Consider both pairwise interaction analysis and higher-order niche definitions.
- Explicitly compare contact-scale and larger neighborhood-scale analyses.

### Clinical association

- Patient-level aggregation is critical when multiple ROIs exist.
- Kaplan-Meier plus Cox modeling recurred across several outcome-linked papers.
- Correct for multiple testing by FDR at minimum.

## 13. Checklist of aspects to investigate in a new dataset

- Cohort:
  - number of patients
  - number of analyzable ROIs after QC
  - outcome maturity
  - clinical covariates
- Assay:
  - platform
  - panel size
  - marker coverage
  - image resolution
  - batch structure
- Segmentation:
  - tool
  - nucleus versus membrane assumptions
  - expansion rules
  - failure modes
- Feature extraction:
  - intensity summarization
  - area normalization
  - outlier handling
  - batch correction
- Cell typing:
  - clustering algorithm
  - annotation strategy
  - rare-state detectability
  - cross-sample reproducibility
- Spatial modeling:
  - contact graph
  - radius-based neighborhoods
  - border distance
  - niche clustering
  - permutation framework
- Validation:
  - independent cohort
  - orthogonal modality
  - pathology concordance
  - parameter sensitivity
- Reporting:
  - code links
  - data links
  - ROI metadata
  - segmentation masks
  - clinical metadata availability

## 14. Gaps worth refining together

This draft synthesis still needs decisions from us on:

- what should count as the default spatial unit for prostate datasets
- which normalization strategy should be considered baseline
- which stromal / immune / epithelial marker sets are minimally required
- which niche-definition approach is most appropriate for prostate tissue:
  - contact graph
  - local composition
  - segmentation-free embedding
- which validation standards should be mandatory versus optional
- how to convert these cross-paper observations into a decision tree for new datasets
