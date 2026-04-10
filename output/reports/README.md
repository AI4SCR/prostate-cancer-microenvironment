# Paper Reports

This folder contains one report per paper from `resources/`, normalized to the same sections:

- data and metadata
- biological key findings
- methods used to obtain the findings
- code and data links
- uncertainties / missing metadata

## Summary table

| Paper | Disease | Technology | Patients | Samples / regions | Cells / objects | Report |
| --- | --- | --- | ---: | --- | --- | --- |
| Keren et al. 2018 | Triple-negative breast cancer | MIBI-TOF | 41 | 42 cores, 1,763 images | 3,000-10,000 cells per patient reported; total not stated | [keren_2018_tnbc_mibi.md](./keren_2018_tnbc_mibi.md) |
| Keren et al. 2019 | Triple-negative breast cancer | MIBI-TOF | 1 | 8 high-resolution regions from 520 tiled FOVs | 11,120 | [keren_2019_mibi_tof_platform.md](./keren_2019_mibi_tof_platform.md) |
| Danenberg et al. 2022 | Breast cancer | IMC | 693 | 749 images from 794 tissue spots | Total not stated | [danenberg_2022_breast_tme_structures.md](./danenberg_2022_breast_tme_structures.md) |
| Jackson et al. 2020 | Breast cancer | IMC | 352 | 720 images | 1,267,078 combined across cohorts | [jackson_2020_single_cell_pathology_breast_cancer.md](./jackson_2020_single_cell_pathology_breast_cancer.md) |
| Bonollo et al. 2020 | Prostate cancer | Review article | Not applicable | Not applicable | Not applicable | [bonollo_2020_cafs_prostate_review.md](./bonollo_2020_cafs_prostate_review.md) |
| Cords et al. 2023 | Breast cancer; cross-cancer validation | scRNA-seq and IMC | 14 discovery; 12 IMC | 7-13 IMC regions per patient | 16,704 stromal cells in discovery; 553,121 IMC objects | [cords_2023_caf_classification.md](./cords_2023_caf_classification.md) |
| Amitay et al. 2024 | Melanoma sentinel lymph nodes | MIBI-TOF and CosMx | 69 | 177 MIBI images; 77 CosMx images | 1,652,274 MIBI cells; 848,727 CosMx cells | [amitay_2024_melanoma_sln.md](./amitay_2024_melanoma_sln.md) |
| Tan et al. 2025 | Lung adenocarcinoma | IMC with GeoMx/scRNA-seq validation | 416 discovery; 143 validation | 76,963 LTME tiles; 96 GeoMx cores | LTME objects reported; total human cells not stated | [tan_2025_canvas_lung_heterogeneity.md](./tan_2025_canvas_lung_heterogeneity.md) |
| Martinelli et al. draft | Primary prostate cancer | IMC | 190 final evaluable; 196 initially reported | 459 tumor ROIs; 523 acquired; 541 in abstract | 2,191,967 | [martinelli_2026_pca_draft.md](./martinelli_2026_pca_draft.md) |

## Notes

- Counts were taken from the paper text when explicitly stated.
- When a paper is a review or does not report a metric directly, the report says so instead of inferring it.
- Some papers report modality-specific counts rather than a single deduplicated total.
- The PCA draft report preserves draft-specific inconsistencies and placeholder sections because they are part of the current manuscript state.
