## Task

Spend up to 24 hours exploring and analyzing the prostate cancer microenvironment (PCa) dataset.

Your goal is to understand what has already been done in this project, identify promising biological questions, and follow the most interesting signals in the data with curiosity and initiative. Prioritize genuinely interesting biological observations, hypotheses, and directions for follow-up over reproducing the full historical pipeline.

## Hard Constraints

You may only modify files inside:

- `$HOME/projects/prostate-cancer-microenvironment/scripts/99-exploratory-analysis`
- `$HOME/projects/prostate-cancer-microenvironment/src`

You are strictly forbidden from modifying any file outside those two directories.

## What To Read First

Before doing analysis, read:

- `$HOME/projects/prostate-cancer-microenvironment/CLAUDE.md`
- `$HOME/projects/prostate-cancer-microenvironment/README.md`
- `$HOME/projects/prostate-cancer-microenvironment/objectives.md`
- the generated reports in this repository

Use these to understand the local instructions, project goals, available tooling, prior outputs, and biological framing.

Then review what analyses have already been done in this repository and, where useful, compare that with relevant prior publications before deciding what to explore more deeply.

## Project Context

- `scripts/99-exploratory-analysis/main.py` is a very short starter example showing how to load the PCa data and begin working with ATHENA-related analysis.
- Treat `main.py` as orientation only, not as a complete pipeline.
- The dataset is accessed through `ai4bmr_datasets.PCa`, typically using `BASE_DIR` from the local `.env`.
- Some packages are installed manually as described in the README.
- If you need additional packages, you may install them, but only using `uv`.

## Suggested Approach

1. Familiarize yourself with the reports and previous outputs to understand what biology and analyses are already covered.
2. Inspect the available data structure, metadata, marker panel, image-level content, and cell-level content.
3. Identify underexplored or promising questions, especially around tumor microenvironment organization, stromal-immune-epithelial structure, heterogeneity, niche-like patterns, and clinically relevant variation across patients or samples.
4. Follow the most interesting directions iteratively. If you find a potentially meaningful signal, test it from multiple angles and try to understand whether it is likely biological, technical, cohort-driven, or confounded.
5. Leave behind readable code, notes, and outputs inside the allowed directories so another researcher can understand what you explored, what you found, and what seems worth pursuing next.

## Quality Bar

- Prefer a small number of thoughtful, well-supported findings over many shallow plots.
- Be explicit about uncertainty and limitations.
- Distinguish clearly between exploratory observations and stronger conclusions.
- When possible, relate findings back to prior analyses in this repository so the novelty is clear.
