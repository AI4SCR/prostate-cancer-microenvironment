# Prostate Cancer Spatial Omics Project

## Overview
This project focuses on constructing a spatially resolved single-cell atlas of primary prostate cancer (PCa) by integrating high-dimensional imaging data with computational analysis and machine learning.

The primary goal is to characterize the tumor microenvironment (TME), quantify inter- and intra-patient heterogeneity, and identify clinically relevant spatial patterns that extend beyond traditional pathology metrics such as Gleason grade.

## Working Instructions

All new exploratory work should be carried out in `scripts/99-exploratory-analysis`.

The file `scripts/99-exploratory-analysis/main.py` is intentionally a very short example that shows how to access the data and start working with ATHENA-related analysis in this project. Treat it as a starting point, not as a complete pipeline.

## Project Objectives

- Build a high-resolution spatial atlas of prostate cancer tissue
- Quantify cell type composition and spatial organization
- Study inter- and intra-patient heterogeneity
- Identify stromal–immune–epithelial niches with prognostic relevance
- Derive image- and sample-level scores that stratify patients by outcome
- Develop scalable computational pipelines for spatial omics
- Move beyond handcrafted features toward representation learning-based approaches

## Data Access

The dataset can be loaded as follows:

```python
ds = PCa(
    base_dir=Path(os.environ["BASE_DIR"]),
    image_version="filtered",
    mask_version="annotated",
    load_metadata=True,
    load_intensity=True,
    align=True,
)
ds.setup(engine="pyarrow")
