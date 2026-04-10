# Spatial Profiling of the Tumor Microenvironment in Prostate Cancer

- Source PDF: `resources/pca-paper.pdf`
- Citation: Martinelli et al., draft manuscript, PDF created April 10, 2026
- Venue / journal: Not stated in the draft
- DOI: Not found

## Data and metadata

- Disease / indication: Primary prostate cancer after radical prostatectomy, with clinical endpoints including biochemical recurrence, clinical progression, disease progression, and overall survival.
- Tissue / specimen type: FFPE prostate TMAs from radical-prostatectomy specimens, with up to four cores per patient sampled from high- and low-Gleason areas.
- Spatial or single-cell technology: Imaging mass cytometry (IMC) with spatial single-cell analysis.
- Marker panel / assay details: Prostate-tailored IMC panel described as 34-plex in the abstract and 35-plex in the methods; one marker, FAP, was later excluded after QC, leaving an effective 34-marker analysis. The panel spans epithelial, stromal / CAF, smooth-muscle, endothelial, immune, and functional markers including ERG, p53, beta-catenin, and pYAP1.
- Number of patients / samples / regions / images: The abstract reports 196 patients and 541 ROIs. Methods and results report 196 initial patients, 190 final evaluable patients, 523 acquired ROIs, and 459 tumor-containing ROIs after filtering.
- Number of cells / objects analyzed: 2,191,967 single-cell profiles.
- Cohorts / validation sets: No external validation cohort is described; validation is internal and includes QC filtering, repeated niche-stability runs, manual niche review, and pathological review of serial H&E sections.

## Biological key findings

- The draft defines a spatial single-cell atlas of primary prostate cancer with 34 cell types organized into 18 recurrent spatial niches.
- These niches group into epithelial-dominant, CAF-dominant, and immune-enriched supergroups.
- An ERG-positive / p53-positive luminal epithelial state is associated with worse overall and progression-free survival.
- Stromal heterogeneity separates ECM-remodeling CD105high CAF1 / myCAF-like states from contractile SMC-like CAF2 states, with additional CES1, EGR1, CD146, and AR-related distinctions.
- A small periglandular CD105high CAF1 niche shows the strongest adverse clinical association and high stromal-immune connectivity.
- Immune niches, including a TLS-like B-cell / T-cell aggregate, track with histological inflammation.
- A composite immune niche burden stratifies worse patient survival, suggesting cumulative immune architecture matters more than any single immune population.
- The draft argues that spatial stromal-immune-epithelial organization captures clinically relevant biology beyond Gleason grade alone.

## Methods used to obtain the findings

- Experimental methods: FFPE TMA sections from radical-prostatectomy samples were stained with an IMC panel; serial H&E sections were used for histopathology annotation; antibodies were pre-validated by immunofluorescence on prostate tissue.
- Preprocessing / segmentation / annotation: Steinbock for hot-pixel removal, CATALYST for spillover assessment, Steinbock plus DeepCell for nuclear segmentation, and per-cell intensities clipped at the 99.9th percentile, arcsinh-transformed, and min-max normalized. Cell phenotyping used iterative k-NN graph construction with Parametric UMAP and Leiden clustering followed by hierarchical reclustering.
- Statistical / computational analysis: CLR-transformed cell-type proportions, Jensen-Shannon ROI distances, hierarchical clustering, univariate Cox proportional hazards models, FDR correction, Kaplan-Meier / log-rank tests, and patient-level max-pooling across ROIs.
- Spatial analysis: Radius-based cell-cell graphs with a 32-pixel neighborhood, k-means clustering on local composition vectors, 50-run stability checks using adjusted Rand index, manual merging into 18 niches, z-score niche enrichment, Spearman niche co-occurrence, and radius-graph interaction profiling.
- Validation: QC to remove artifacts and segmentation errors, repeated k-means stability checks, manual biological validation of niches, and comparison of spatial scores against H&E-derived stromogenic and inflammatory annotations.

## Code and data availability

- GitHub / code: Not reported in the draft text reviewed. The code-availability section still contains placeholder text: `Adriano: github link`.
- Data: Not reported in the draft text reviewed. A data-availability section exists but appears unfinished.

## Draft-specific caveats / missing metadata

- The manuscript contains multiple placeholder citations marked `REF`.
- Cohort counts are inconsistent across sections: abstract says 196 patients / 541 ROIs, while methods and results describe 190 final evaluable patients, 523 acquired ROIs, and 459 tumor-containing ROIs.
- Panel size is also inconsistent: 34-plex in the abstract versus 35-plex in methods, with FAP later excluded after QC.
- Data-availability and code-availability sections are incomplete.
- TLS-like structures are inferred from spatial organization only; the draft notes that dedicated TLS maturation / localization markers were not included.
- Stromogenic status is based on non-standardized histopathologic criteria.
- The cohort has relatively few clinical events, which the draft notes as a limitation for survival power.

```yaml
machine_readable:
  title: "Spatial Profiling of the Tumor Microenvironment in Prostate Cancer"
  year: 2026
  disease: "primary prostate cancer"
  technology: "imaging mass cytometry (IMC)"
  n_patients: 190
  n_samples: 459
  n_cells: 2191967
  code_links: []
  data_links: []
```
