# Breast tumor microenvironment structures are associated with genomic features and clinical outcome

- Source PDF: `resources/s41588-022-01041-y.pdf`
- Citation: Danenberg et al., 2022, *Nature Genetics*
- DOI: `10.1038/s41588-022-01041-y`

## Data and metadata

- Disease / indication: Breast cancer, specifically treatment-naive primary tumors from the METABRIC cohort.
- Tissue / specimen type: FFPE breast tumor TMAs, mostly 0.6 mm spots with some 1 mm spots; 4 micron sections were profiled.
- Spatial or single-cell technology: Imaging mass cytometry (IMC) with multitiered spatial network analysis.
- Marker panel / assay details: 37-dimensional IMC images; the study defined 32 final cell phenotypes, split into 16 epithelial and 16 TME phenotypes.
- Number of patients / samples / regions / images: 794 tissue spots from 718 patients were processed initially; after excluding 31 normal and 14 in situ-only spots, 749 images from 693 patients remained. Of these tumors, 635 were represented by one spot, 55 by two spots, and 3 by three spots.
- Number of cells / objects analyzed: A single cohort-wide total cell count was not explicitly reported in the paper text reviewed.
- Cohorts / validation sets: Discovery clustering used 458 tumors from one METABRIC center, and validation used 181 tumors from the second center through a trained random-forest classifier.

## Biological key findings

- The study identified 10 recurrent multicellular TME structures spanning quiescent and active stroma, vascularized stroma, active immune-response states, TLS-like structures, and a suppressed-expansion state.
- Suppressed-expansion and TLS-like structures were the largest and most diverse structures, with some containing more than 200 cells.
- B cells became more dominant as local connectivity increased in several structures, while T-cell proportions remained comparatively stable, implying distinct organization principles for lymphoid compartments.
- TME structure was strongly associated with genomic subtype, and network-derived features were especially predictive of IntClust 4- tumors with out-of-sample AUC 0.91.
- Suppressed-expansion structures were enriched for BRCA1 and CASP8 mutations, linking immune suppression and DNA-repair-related genomic states.
- Additional genomic associations included BRCA2 with active immune-response and vascular / stromal structures, CD274 gains with granulocyte-enriched structures, B2M loss with TLS-like structures, and CDH1 mutations with active stromal / vascular structures.
- In ER-positive disease, granulocyte-enriched, APC-enriched, and suppressed-expansion structures were associated with poor survival, while vascular stroma was associated with better outcome.

## Methods used to obtain the findings

- Experimental methods: FFPE METABRIC TMAs stained with a 37-marker IMC panel and ablated on an imaging mass cytometer; linked genomic and clinical annotations were incorporated from METABRIC.
- Preprocessing / segmentation / annotation: Ilastik pixel classification, CellProfiler segmentation, CATALYST spillover compensation, Gaussian-mixture modeling of pan-cytokeratin for epithelial separation, and SOM plus PhenoGraph clustering with manual phenotype merging.
- Statistical / computational analysis: Shannon diversity, generalized linear models, Ward hierarchical clustering, consensus clustering to define 10 TME structures, random-forest validation, regularized logistic regression for subtype prediction, Cox models, log-rank tests, and Benjamini-Hochberg correction.
- Spatial analysis: Separate epithelial and non-epithelial cell-contact graphs, community detection on spatial graphs, connectivity-profile summaries, and network metrics including diameter, density, transitivity, and assortativity.
- Validation: Independent-center validation using a trained random forest, out-of-sample subtype prediction across METABRIC centers, ER-stratified survival analysis, and genomic association analysis against public METABRIC resources.

## Code and data availability

- GitHub / code:
  - `https://github.com/BodenmillerGroup/ImcSegmentationPipeline`
  - `https://github.com/CellProfiler`
  - `https://github.com/ilastik`
  - `https://github.com/saeyslab/FlowSOM`
  - `https://github.com/JinmiaoChenLab/Rphenograph`
  - `https://github.com/igraph/rigraph`
  - `https://doi.org/10.5281/zenodo.6036188`
- Data:
  - `https://doi.org/10.5281/zenodo.5850952`
  - `https://www.cbioportal.org/`
  - `EGAS00000000083`
  - `EGAS00001001753`

## Uncertainties / missing metadata

- The paper does not state one total cohort-wide cell count in the main text reviewed.
- The sample count is layered: 794 spots were processed initially, but 749 images from 693 patients remained after filtering.
- The discovery / validation splits are clear for some analyses, but not every downstream analysis uses exactly the same subset.
- Some clinical covariates and genomic details are in linked resources or supplementary material rather than the main text.

```yaml
machine_readable:
  title: "Breast tumor microenvironment structures are associated with genomic features and clinical outcome"
  year: 2022
  disease: "breast cancer"
  technology: "imaging mass cytometry (IMC)"
  n_patients: 693
  n_samples: 749
  n_cells: null
  code_links:
    - "https://github.com/BodenmillerGroup/ImcSegmentationPipeline"
    - "https://github.com/CellProfiler"
    - "https://github.com/ilastik"
    - "https://github.com/saeyslab/FlowSOM"
    - "https://github.com/JinmiaoChenLab/Rphenograph"
    - "https://github.com/igraph/rigraph"
    - "https://doi.org/10.5281/zenodo.6036188"
  data_links:
    - "https://doi.org/10.5281/zenodo.5850952"
    - "https://www.cbioportal.org/"
    - "EGAS00000000083"
    - "EGAS00001001753"
```
