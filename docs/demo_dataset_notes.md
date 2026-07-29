# Note on the MIMIC-IV Demo prototyping pass

**Status:** PhysioNet credentialed access to full MIMIC-IV is still pending manual
review (automated credentialing rejected the application over the academic-email
requirement; appeal sent to credentialing@physionet.org citing the ORCID
alternative permitted by PhysioNet's own FAQ).

Rather than wait idle, `scripts/build_cohort_demo.py` builds and validates the
Phase 2 cohort-extraction logic against the **MIMIC-IV Clinical Database Demo**
— an openly available, no-credentialing-required 100-patient subset of MIMIC-IV,
published by the same team (Johnson et al., 2023, DOI: 10.13026/dp1f-ex47).
Same schema, same table structure as the full database.

## What this validates

- The full join chain across `patients`, `admissions`, `icustays`, `diagnoses_icd`,
  `labevents` (lactate), `chartevents` (MAP), and `inputevents` (vasopressors)
  runs correctly end-to-end.
- T0 identification logic (first lactate > 4.0 mmol/L with a concurrent MAP >
  65 mmHg, per PROTOCOL.md Section 3) executes without errors.
- Treatment-group assignment (vasopressor within the 3-hour grace period vs.
  standard care) and the 28-day mortality outcome calculation both run correctly.

## What's intentionally simplified for this prototype (must change before final analysis)

1. **Sepsis-3 proxy:** the demo script flags "sepsis" using a static list of
   ICD-9/ICD-10 sepsis diagnosis codes on the hospital admission, rather than
   the formal Sepsis-3 (Seymour et al.) algorithm — suspected infection
   (antibiotic + culture timing) plus SOFA ≥ 2. The full pipeline should use
   the official `sepsis3.sql` concept table from the
   [MIT-LCP/mimic-code](https://github.com/MIT-LCP/mimic-code) repository.
2. **No SOFA / Charlson Comorbidity Index computed yet** — PROTOCOL.md's
   propensity score model (Section 4.1) requires both; these should also be
   pulled from mimic-code's concept tables against the full database, not
   reimplemented from scratch.
3. **Exclusions not yet implemented:** active hemorrhage, DNR within 24h of
   admission (PROTOCOL.md Section 3, Component 1). Deferred until full-data
   pass since they weren't needed to validate the join logic.
4. **Sample size is not meaningful.** With only 100 total patients in the
   demo, expect a single-digit (or zero) eligible cohort — this run found
   1 qualifying stay. This is expected and does not indicate a bug; it
   reflects the demo's small, random subject sample. No inference should be
   drawn from demo-run counts or outcomes.

## Migration path to full data

Once credentialed access is approved, port `build_cohort_demo.py`'s query
logic to the BigQuery-hosted full MIMIC-IV database (per PROTOCOL.md Section
5), swap in the mimic-code concept tables for Sepsis-3/SOFA/Charlson, add the
missing exclusion criteria, and re-run. The join structure and T0/treatment
assignment logic validated here should not need to change.
