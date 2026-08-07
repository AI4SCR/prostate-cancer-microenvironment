# prostate-cancer-microenvironment

Code companion to *"Spatial single-cell proteomics defines multicellular niches
in the primary prostate cancer microenvironment"* (bioRxiv
[10.64898/2026.04.30.721907v1](https://doi.org/10.64898/2026.04.30.721907)).
This repo reproduces the paper's figures from already-processed single-cell
data — one script per figure/panel, each a verbatim (or explicitly disclosed)
port of the original analysis code.

**Start here:**
- [`figures.md`](figures.md) — the authoritative figure → script mapping,
  with per-panel validation status and known issues. Check this before
  running anything for a specific figure.
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) — the full pipeline from raw
  acquisitions to labeled cells, environment troubleshooting, and every known
  discrepancy against the published figures.
- [`CLAUDE.md`](CLAUDE.md) — conventions for anyone (human or agent) adding
  or modifying scripts here (verbatim-port discipline, no new-data-computing
  rule, etc.).

> **Note on data provenance**: this repo currently reads from data staged
> locally on this machine (`BASE_DIR`, `LEGACY_DATA_DIR` below) rather than
> downloading a single self-contained archive. A consolidated Zenodo release
> that makes the whole pipeline runnable end-to-end from one download
> (`zenodo download → uv sync → run scripts`) is planned but not yet
> published — for now, ask a project maintainer for access to the staged
> data directories if you don't already have them.

## Project structure

```text
prostate-cancer-microenvironment/
├── scripts/
│   ├── data/              # export.py (raw → EXPORT_DIR tables) + UMAP reducer ports
│   ├── figures/           # one script per figure/panel — see figures.md
│   └── port/              # port_umap_reducer.py, its own pinned pixi env
├── src/prostate_cancer/    # shared utilities (colormaps, normalization, plotting helpers)
├── resources/              # colormaps.yaml, metalabels.yaml, non_marker_channels.txt
├── data/                   # generated + staged legacy tables (gitignored, machine-specific)
├── output/figures/         # generated figure PDFs/PNGs (gitignored)
├── figures.md               # figure → script mapping + validation status
├── REPRODUCIBILITY.md       # full pipeline, environment notes, known discrepancies
└── CLAUDE.md                 # contribution conventions
```

## Setup

### 1. Configure `.env`

```bash
cp .env.example .env
```

Fill in the machine-specific paths (see comments in `.env.example` for what
each one is for and why it's separate from the others):

- `BASE_DIR` — root of the staged `ai4bmr-datasets` `PCa` layout (read-only,
  reproduces the published Zenodo dataset `10.5281/zenodo.19552665`).
- `EXPORT_DIR` — where generated intermediate tables go (defaults to
  `data/` inside this repo).
- `OUTPUT_FIGURES_DIR` — where generated figure files go (defaults to
  `output/figures/` inside this repo).
- `LEGACY_DATA_DIR` — consolidated, in-repo copy of pre-migration legacy data
  that has no reproducing script here yet (niche-neighborhood graph data,
  precomputed niche annotation/composition tables). See `REPRODUCIBILITY.md`
  for what's in it and why.

### 2. Python environment

```bash
uv sync
```

This installs everything in `pyproject.toml`, including `ai4bmr-datasets`
(pinned to its `pca` branch — the `PCa` dataset class this repo depends on
isn't on `main` yet).

Run any Python figure script with:

```bash
uv run python scripts/figures/figure4_patient_clustering.py
```

### 3. R environment

R scripts load the same `.env` via `dotenv::load_dot_env()`. Required
packages: `dotenv`, `arrow`, `tidyverse`, `survival`, `survminer`,
`gtsummary`, `compositions`, `coxme`, `ComplexHeatmap`, `circlize`, `yaml`,
`entropy`. No `renv.lock` is maintained — versions aren't pinned on the R
side. See `REPRODUCIBILITY.md`'s Setup section for a `survminer` install
workaround needed on some R 4.4.1 builds.

On this cluster, load R via the module system before running any R script:

```bash
module load r-light/4.4.1
Rscript scripts/figures/figure2_cell_type_heatmap.R
```

### 4. UMAP reducer porting (one-time, only if regenerating UMAP-based figures)

Figures 2b, 3a, S2, and S7 read pre-computed UMAP embeddings rather than
re-fitting UMAP (the legacy pipeline never seeded `UMAP.fit()`, so a fresh
fit is not reproducible against the published panels). These embeddings are
ported once from legacy `reducer.pkl` files using a separately pinned
environment (older `numba`/`umap-learn` stack required to unpickle them):

```bash
cd scripts/port
pixi run python port_umap_reducer.py <path/to/reducer.pkl> <out_path.parquet>
```

See each figure's own `scripts/data/figureN_*_embedding.py` docstring for
the exact `reducer.pkl` paths and output locations — the ported parquet
files are cached under `data/figures/`, so this only needs to run once per
config. `figures.md` lists which figures need this.

## Running the figures

1. **Generate the base tables** (`metadata.parquet`, `clinical.parquet`,
   `intensity.parquet`, `intensity_normalized.parquet`), once, before
   anything else:

   ```bash
   uv run python scripts/data/export.py
   ```

2. **Look up the script for the figure/panel you want** in
   [`figures.md`](figures.md)'s Overview table — it lists the exact script,
   its inputs, its output path, and whether it's already been validated
   against the published panel.

3. **Run it.** Python scripts:

   ```bash
   uv run python scripts/figures/figure6c_niche_composition_filtered.py
   ```

   R scripts:

   ```bash
   module load r-light/4.4.1
   Rscript scripts/figures/figure6_niche_abundance_heatmap.R
   ```

   Every script writes to `$OUTPUT_FIGURES_DIR/figureN/...` — see
   `figures.md`'s `Output` column for the exact filename(s).

4. **To run every figure**, loop over `scripts/figures/*.py` /
   `scripts/figures/*.R` — there's no single "run everything" entrypoint
   yet, since some scripts (UMAPs, KM survival curves over the full cohort)
   are heavy enough to warrant a SLURM job rather than an interactive login
   shell (see `AGENTS.md` if you're running on the HPC cluster).

## Raw data acquisition (advanced, rarely needed)

`BASE_DIR` above is expected to already contain segmented, feature-extracted
cells (the published Zenodo dataset snapshot). The steps below — needed only
if you're regenerating that snapshot from raw `.mcd` acquisitions — are
summarized here; the full pipeline (clustering, label transfer, bootstrap
loop) is documented in `REPRODUCIBILITY.md`'s "Pipeline: raw acquisitions →
labeled cells" section.

### Steinbock (segmentation + feature extraction)

Install [Docker Desktop](https://docs.docker.com/desktop/install/mac-install/)
and make sure the daemon is running, then add this to your shell rc file:

```bash
steinbock0.16.1(){
docker run \
-v "$(pwd)":/data \
ghcr.io/bodenmillergroup/steinbock:0.16.1 \
"$@"
}
```

```bash
mkdir -p steinbock/raw && cd steinbock
find <PATH_TO_RAW_DATA> -name "*.mcd" -exec cp {} raw/ \;
ln -s <PATH_TO_PROJECT>/03_spatial/03_spatial/process-steinbock-panel.py .
steinbock0.16.1 preprocess imc panel

conda activate PCa
python process-steinbock-panel.py panel.csv

steinbock0.16.1 preprocess imc images --hpf 50 --imgout img
steinbock0.16.1 segment deepcell --minmax --type nuclear && \
steinbock0.16.1 measure intensities && \
steinbock0.16.1 measure regionprops && \
steinbock0.16.1 measure neighbors --type centroids --kmax 5
```

### R install note (Linux, non-Homebrew)

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.18")
```

## Troubleshooting

- Symlinks: make sure the target of any symlink is available inside the
  Docker image (i.e. actually mounted), not just on the host.
- R/UMAP/`ai4bmr-datasets` environment issues: see `REPRODUCIBILITY.md`'s
  Setup section first — several non-obvious version-pinning issues (e.g.
  `survminer`'s `Deriv` dependency) are already documented there with fixes.
