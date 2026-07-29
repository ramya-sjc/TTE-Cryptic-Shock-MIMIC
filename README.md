# Target Trial Emulation: Vasopressor Timing in Cryptic Shock (MIMIC-IV)

**Status: In Progress - Phase 1 (Protocol complete, PhysioNet credentialing underway)**

A real-world evidence (RWE) study replicating and extending a February 2026 published target trial emulation, using ICU electronic health record data (MIMIC-IV) to estimate the causal effect of early vs. delayed vasopressor initiation on 28-day mortality in patients with "Cryptic Shock" (severe tissue hypoperfusion despite normal blood pressure).

Built as a portfolio project applying pharmacoepidemiology and causal inference methods to real-world clinical data - the same methodological framework used by real-world evidence teams in pharmaceutical and healthcare settings.

---

## Why this project

Naive comparisons of "early vs. late" treatment in observational EHR data are biased by two well-known problems: confounding by indication (sicker patients get treated sooner) and immortal time bias (patients must survive long enough to be classified into the "late" group). Target trial emulation (Hernán & Robins, 2016) addresses both by explicitly specifying the hypothetical randomized trial being emulated *before* looking at outcomes, then applying causal inference methods (here, inverse probability of treatment weighting) to approximate it.

This project first **replicates** the design and analysis of a recently published paper (Shuaibu et al., 2026, *Research Square* preprint) as a validation check on my own pipeline, then **extends** it with an original analysis addressing a limitation the original authors flagged themselves. See [`PROTOCOL.md`](./PROTOCOL.md) for the full seven-component target trial specification, statistical analysis plan, and literature review.

---

## Repository structure

```
.
├── README.md                  # this file
├── PROTOCOL.md                 # full target trial protocol (7 components + analysis plan)
├── scripts/
│   └── build_cohort_demo.py    # cohort-extraction prototype, run against the open MIMIC-IV Demo
├── sql/                        # cohort extraction queries (Phase 2, pending full data access)
│   └── (to be added)
├── analysis/                   # R scripts: propensity scores, IPTW, survival analysis (Phase 3-4)
│   └── (to be added)
├── report/                     # final rendered report, tables, figures (Phase 6)
│   └── (to be added)
└── docs/
    ├── literature/              # annotated notes on cited papers
    │   └── (to be added)
    └── demo_dataset_notes.md    # what the demo prototype validates vs. what's still simplified
```

---

## Roadmap

- [x] **Phase 0 — Setup:** PhysioNet account created; CITI "Data or Specimens Only Research" training in progress
- [x] **Phase 1 — Protocol:** Target trial protocol drafted and version-controlled (`PROTOCOL.md`, v1.0, 2026-07-13)
- [x] **Phase 2 (prototype) — Cohort extraction logic validated** against the openly-available MIMIC-IV Clinical Database Demo (100 patients, no credentialing required) — see `scripts/build_cohort_demo.py` and `docs/demo_dataset_notes.md`
- [ ] **Phase 2 (full) — Cohort extraction:** same logic re-run against full MIMIC-IV (BigQuery), pending credentialed access approval
- [ ] **Phase 3 — Trial emulation design:** Strategy assignment, immortal-time-bias handling (cloning/censoring/weighting)
- [ ] **Phase 4 — Statistical analysis:** IPTW, weighted Cox model, balance diagnostics
- [ ] **Phase 5 — Sensitivity & extension analysis:** Robustness checks + original negative-control (or alternative) analysis
- [ ] **Phase 6 — Report & write-up:** Final reproducible report (Quarto/R Markdown → PDF/HTML)

Current blocker: PhysioNet credentialing application still awaiting manual review (appealed to credentialing@physionet.org). Cohort-extraction pipeline is being built and validated now against the open MIMIC-IV Demo so no time is lost waiting.

---

## Data access

This project uses **MIMIC-IV v3.1**, a restricted-access, de-identified ICU database hosted on [PhysioNet](https://physionet.org), requiring individual credentialed access (CITI training + signed Data Use Agreement). Per PhysioNet's data use terms, raw data is never included in this repository - only code, protocol documents, and de-identified/aggregated outputs (tables, figures, summary statistics) will be published once analysis is complete.

---

## Methods summary

- **Population:** Adult ICU patients meeting Sepsis-3 criteria with "Cryptic Shock" (MAP > 65 mmHg + lactate > 4.0 mmol/L)
- **Exposure:** Vasopressor initiation within a 3-hour grace period from time-zero, vs. standard fluid-only care
- **Outcome:** 28-day all-cause mortality (primary); ICU length of stay and a negative control outcome (secondary/extension)
- **Confounding adjustment:** Inverse Probability of Treatment Weighting (IPTW) via propensity score (age, SOFA score, Charlson Comorbidity Index, lactate, MAP, pre-T0 fluid volume)
- **Analysis:** Weighted Kaplan-Meier curves, IPTW-weighted Cox proportional hazards model, standardized mean difference balance diagnostics

Full detail in [`PROTOCOL.md`](./PROTOCOL.md).

---

## Key references

- Shuaibu II, Radouani S, Hussain Y. Early Vasopressor Initiation versus Standard Fluid Resuscitation in Normotensive Cryptic Shock: A Target Trial Emulation. *Research Square* [Preprint], 2026. DOI: [10.21203/rs.3.rs-8747647/v1](https://doi.org/10.21203/rs.3.rs-8747647/v1)
- Hernán MA, Robins JM. Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available. *Am J Epidemiol.* 2016;183(8):758-764.
- Johnson AEW, Bulgarelli L, Pollard TJ, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data.* 2023;10(1):1.

---

## Author

Ramya Sri Jayashanker Chithra — MSc Health Data Science
Portfolio project for real-world evidence / pharmacoepidemiology roles.

## License

Code and documentation in this repository are shared under the MIT License. This repository does not and will not contain any raw MIMIC-IV patient data, per PhysioNet's Data Use Agreement.
