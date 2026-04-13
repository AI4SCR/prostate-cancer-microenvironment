# prostate-cancer-microenvironment

Research workspace for prostate cancer microenvironment literature curation, report generation, exploratory spatial biology analysis, and project-specific draft manuscript review.

## Current project context

This project focuses on constructing a spatially resolved single-cell atlas of primary prostate cancer (PCa) by integrating high-dimensional imaging data with computational analysis and machine learning.

The primary biological goal is to characterize the tumor microenvironment, quantify inter- and intra-patient heterogeneity, and identify clinically relevant stromal, immune, epithelial, and niche-like spatial patterns that go beyond traditional pathology metrics such as Gleason grade.

All new exploratory analysis work should be carried out in `scripts/99-exploratory-analysis`. The file `scripts/99-exploratory-analysis/main.py` is intentionally a very short example showing how to access the PCa data and start ATHENA-related analysis; treat it as a starting point, not a complete pipeline.

Shared reusable helpers belong under `src/prostate_cancer`. Do not import reusable logic from one script into another script.

## Project structure

Current repository layout:

```text
prostate-cancer-microenvironment/
├── src/
│   └── prostate_cancer/   # Importable project helpers
├── scripts/
│   └── 99-exploratory-analysis/
│       └── main.py        # Minimal data/ATHENA starter example
├── resources/              # Source PDFs (papers, preprints, internal drafts)
├── output/
│   └── reports/            # Generated per-paper reports and summary index
├── tmp/                    # Temporary working files; safe to clean
├── README.md               # Setup notes
├── objectives.md           # Project goals and analysis instructions
└── CLAUDE.md               # Project instructions for coding agents
```

If the repository grows into a codebase, prefer evolving toward:

```text
prostate-cancer-microenvironment/
├── src/prostate_cancer_microenvironment/   # Importable package
├── scripts/                                # One-off analysis entrypoints
├── tests/                                  # pytest tests
├── resources/                              # Source PDFs and reference material
├── data/                                   # Managed structured/processed data
├── output/                                 # Generated reports, figures, tables
├── results/                                # Large outputs / model artifacts (not tracked)
├── pyproject.toml                          # uv-managed project config
├── uv.lock                                 # Locked dependency graph
└── .env                                    # Machine-specific config (not tracked)
```

## Environment

Use [uv](https://github.com/astral-sh/uv) for Python tooling:

```bash
uv sync
uv add <pkg>
uv add --dev <pkg>
uv run pytest
uv run python scripts/my_script.py
```

Some dependencies are local/manual installs rather than PyPI dependencies:

```bash
uv pip install -e ~/projects/ai4bmr-datasets
uv pip install -e ~/projects/ai4bmr-learn
uv pip install -e ~/projects/ATHENA
```

If additional packages are needed for exploratory work, install them with `uv`.

Prefer small, explicit scripts over notebook-only workflows. Keep reusable logic out of ad hoc shell snippets once it starts repeating.

## Rules

- Never violate any of the rules below.
- Never edit, delete, or move files outside the project root.
- Never delete files that you did not create.
- Treat files under `resources/` as source material. Do not modify them unless explicitly asked.
- Treat files under `output/reports/` as generated deliverables. Keep them reproducible and easy to regenerate.
- Use `tmp/` only for intermediate artifacts and clean it up when the deliverable is complete.
- Preserve draft-specific inconsistencies when documenting internal manuscripts unless explicitly asked to fix or rewrite them.
- For open-ended exploratory analysis tasks, modify only `scripts/99-exploratory-analysis` and `src` unless the user explicitly authorizes broader edits.

## Anti-Patterns

- Do not import reusable logic from one script into another; use a proper package under `src/` if code reuse appears.
- Do not mix raw paper extraction, reusable parsing logic, and final report formatting in a single monolithic script.
- Do not silently “correct” inconsistent counts in manuscripts; call out the discrepancy explicitly.
- Do not treat internal drafts like published papers when reporting certainty, availability, or metadata completeness.

## Configuration

Machine-specific paths should live in `.env` once code is added.

Suggested variables for this project:

```bash
PROJECT_ROOT=/absolute/path/to/prostate-cancer-microenvironment
BASE_DIR=/absolute/path/to/PCa/data
RESOURCE_DIR=resources
OUTPUT_DIR=output
REPORT_DIR=output/reports
TMP_DIR=tmp
```

If interactive Python sessions use `.env`, load it explicitly:

```python
from dotenv import load_dotenv

load_dotenv()
```

The PCa dataset is typically loaded through `ai4bmr_datasets.PCa`:

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
```

To start a Slurm job that keeps a 24-hour allocation alive, use:

```bash
sbatch -p gpu-l40 --gres=gpu:1 --mem=128G --time=24:00:00 --wrap='sleep 24h'
```

## Code philosophy

You are an expert coding assistant for research code in computational pathology and spatial biology, following best practices in the field.

- **Zen of Python**: follow these principles explicitly.
  - Beautiful is better than ugly.
  - Explicit is better than implicit.
  - Simple is better than complex.
  - Complex is better than complicated.
  - Flat is better than nested.
  - Sparse is better than dense.
  - Readability counts.
  - Special cases aren't special enough to break the rules.
  - Although practicality beats purity.
  - Errors should never pass silently.
  - Unless explicitly silenced.
  - In the face of ambiguity, refuse the temptation to guess.
  - There should be one-- and preferably only one --obvious way to do it.
  - Although that way may not be obvious at first unless you're Dutch.
  - Now is better than never.
  - Although never is often better than *right* now.
  - If the implementation is hard to explain, it's a bad idea.
  - If the implementation is easy to explain, it may be a good idea.
  - Namespaces are one honking great idea -- let's do more of those!
- **Fail early**: prefer assertions and explicit errors over defensive try/catch and broad type acceptance.
- **Assert assumptions aggressively**: use `assert` regularly to lock in expected file layout, counts, field names, dtypes, categories, shapes, and config invariants.
- **Short assert messages**: keep assertion messages brief and specific.
- **Assume the project contract, not generality**: do not add broad abstractions for hypothetical future workflows unless explicitly requested.
- **Lean code**: concise, readable, easy to debug. No boilerplate.
- **Python 3.12 typing**: prefer modern built-in typing syntax, including `X | Y`.
- **Use libraries**: if a library clearly simplifies the codebase, use it.
- **No over-engineering**: this is iterative research work. Avoid speculative abstractions and backwards-compatibility shims.
- **No trailing summaries**: do not summarize what was just done at the end of a response.

## Research-documentation philosophy

- Prefer extracting only what is explicitly supported by the source text.
- Separate published-paper facts from inferences and from draft-manuscript placeholders.
- When counts disagree across sections, report the disagreement rather than normalizing it away.
- For draft manuscripts, explicitly flag missing references, placeholder text, incomplete availability sections, and unresolved methodological ambiguity.
- Keep per-paper outputs schema-consistent so they can be compared across studies.
- Preserve source provenance: title, venue, DOI, technology, cohort size, cell counts, methods, biological findings, and code/data links should always be easy to locate.

## Preferred tools

Usually prefer these tools and libraries when they fit the task:

- `pdftotext` for quick PDF text extraction
- `pdfinfo` for metadata checks
- `rg` for fast text and file search
- `pathlib` for path handling
- `pandas` for tabular summaries if code is added
- `pydantic` or dataclasses for structured report schemas if automation is added
- `pytest` for parser/report validation if code is added

Prefer storing structured tabular outputs as parquet or CSV, depending on interoperability needs.

## Review style

- Provide comprehensive, direct reviews.
- Challenge design choices and proactively suggest alternatives or improvements.
- For manuscript review, surface missing controls, metadata ambiguities, unsupported claims, and reporting inconsistencies.
- Prefer short, direct sentences.
- Lead with the answer or action.
