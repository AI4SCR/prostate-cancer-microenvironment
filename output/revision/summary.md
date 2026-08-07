## Issue 2 follow-up: multivariate Cox results

Following up on the reviewer's question about whether niche 6 and niche 9 are really "independently associated" with outcome "beyond Gleason grade" — we ran the multivariate models they asked for: niche 6/9 together with Gleason grade, and niche 6/9 together with their defining cell-type population, for both overall survival and disease progression.

The short version is that the multivariate results don't back up the "independent, beyond Gleason grade" phrasing as it's currently written. Sharing the numbers below so we can figure out together how best to adjust the text for the response letter.

**Regression formulas**

Every model is a standard Cox proportional-hazards fit, `coxph(Surv(time, event) ~ ..., data)`, with:

- **Outcome (time, event) pair** — one of:
  - Overall survival: `time = last_fu`, `event = (os_status == "dead")`
  - Disease progression: `time = disease_progr_time`, `event = (disease_progr == 1)`
- **Niche term** (`niche6_clr` / `niche9_clr`) — the niche's per-patient abundance, CLR-transformed (centered log-ratio, across the full niche composition) so it's a continuous covariate rather than a high/low split
- **Cell-type term** (`erg_p53_clr` / `caf1cd105_clr`) — the defining cell type's per-patient proportion, CLR-transformed the same way
- **Gleason term** (`gs_grp`) — patient-level ISUP Grade Group (1–5), entered as a plain numeric/ordinal covariate

Models fit:

| Comparison | Univariate | Multivariate |
|---|---|---|
| Niche + Gleason | `Surv(time, event) ~ gs_grp` | `Surv(time, event) ~ niche_clr + gs_grp` |
| Niche + cell type | `Surv(time, event) ~ niche_clr` and `Surv(time, event) ~ celltype_clr` | `Surv(time, event) ~ niche_clr + celltype_clr` |

Each of these is fit once for OS and once for disease progression, giving the 4 + 4 = 8 rows in the tables below.

**Niche + Gleason grade** (patient-level ISUP grade group, n=189, 1 patient missing grade group dropped)

| Niche | Outcome | Niche HR (adjusted for Gleason) | Niche p | Gleason HR | Gleason p |
|---|---|---|---|---|---|
| 6 | OS | 1.16 | 0.11 | 1.39 | 0.014 |
| 6 | Progression | 1.10 | 0.24 | 1.21 | 0.052 |
| 9 | OS | 1.07 | 0.54 | 1.49 | 0.0014 |
| 9 | Progression | 1.14 | 0.11 | 1.24 | 0.020 |

Gleason grade holds up (significant or close to it) in all four fits; the niche term doesn't reach significance in any of them once Gleason is in the model. So on this data, we likely can't claim niche 6 or niche 9 adds prognostic information beyond Gleason grade — worth softening or reframing that line in the text.

**Niche + defining cell type** (n=190)

| Niche | Cell type | Outcome | Niche alone | Niche adjusted | Cell type adjusted |
|---|---|---|---|---|---|
| 6 | ERG+p53+ | OS | HR 1.28, p=0.0033 | HR 0.95, p=0.76 | HR 1.41, p=0.049 |
| 6 | ERG+p53+ | Progression | HR 1.16, p=0.055 | HR 0.85, p=0.24 | HR 1.44, p=0.0062 |
| 9 | CAF1(CD105+) | OS | HR 1.10, p=0.43 | HR 1.34, p=0.17 | HR 0.67, p=0.23 |
| 9 | CAF1(CD105+) | Progression | HR 1.15, p=0.097 | HR 1.36, p=0.042 | HR 0.70, p=0.17 |

Niche 6 is the more interesting case here: it's significant on its own for OS (p=0.0033), but once ERG+p53+ content joins the model its HR moves from 1.28 to 0.95 (p=0.76). That reads as niche 6's OS association largely tracking its ERG+p53+ content, rather than the niche label adding something on top — a nice, explainable result, just not the "independent" one. Same direction, a bit weaker, for progression.

Niche 9 is less clean and doesn't cross significance against CAF1(CD105+) for OS in either direction. For progression it's actually the one result going the other way — not significant alone (p=0.097), just over the line once CAF1 is added (p=0.042) — so this pairing may be worth a second look rather than folding it into the same conclusion as the others.

**Where this leaves us:** across the 8 fits, the niche term is significant on its own in only one (niche 6, OS, univariate), and that doesn't hold once we adjust for ERG+p53+. My read is that the cleanest path forward is to either soften "independently associated, beyond Gleason grade" to something univariate-scoped, or explicitly note that the niche association looks like it's substantially explained by Gleason grade and by the niche's component cell types — happy to help draft either version. Full CSVs with all coefficients are in `output/revision/figure6_niche6/` and `output/revision/figure7_niche9/` if it's useful to look at the raw numbers together.
