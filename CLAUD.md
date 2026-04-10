# prostate-cancer-microenvironment

Research workspace for prostate cancer microenvironment literature curation, report generation, and project-specific draft manuscript review.

## Project structure

Current repository layout:

```text
prostate-cancer-microenvironment/
├── resources/              # Source PDFs (papers, preprints, internal drafts)
├── output/
│   └── reports/            # Generated per-paper reports and summary index
├── tmp/                    # Temporary working files; safe to clean
└── CLAUD.md                # Project instructions for coding agents
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

This repository is not yet set up as a uv-managed Python package.

If Python tooling is added, use [uv](https://github.com/astral-sh/uv) and standardize on:

```bash
uv sync
uv add <pkg>
uv add --dev <pkg>
uv run pytest
uv run python scripts/my_script.py
```

Until then:

- Prefer small, explicit scripts over notebook-only workflows.
- Keep reusable logic out of ad hoc shell snippets once it starts repeating.
- Do not assume package metadata or a managed virtualenv already exists.

## Rules

- Never violate any of the rules below.
- Never edit, delete, or move files outside the project root.
- Never delete files that you did not create.
- Treat files under `resources/` as source material. Do not modify them unless explicitly asked.
- Treat files under `output/reports/` as generated deliverables. Keep them reproducible and easy to regenerate.
- Use `tmp/` only for intermediate artifacts and clean it up when the deliverable is complete.
- Preserve draft-specific inconsistencies when documenting internal manuscripts unless explicitly asked to fix or rewrite them.

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
