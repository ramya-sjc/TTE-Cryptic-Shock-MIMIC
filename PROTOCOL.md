# Target Trial Emulation Protocol

## Early Vasopressor Initiation in "Cryptic Shock": A Replication and Extension Using MIMIC-IV

**Status:** Draft v1.0 - Protocol stage (pre-data-access)
**Author:** Ramya Sri Jayashanker Chithra
**Date:** 2026-07-13
**Repository phase:** Phase 1 - Protocol (see README.md roadmap)

---

## 0. Project status note

This protocol is being written **before** data access is granted, deliberately - this mirrors how real RWE teams work: the analysis plan is locked in before touching the data, so that decisions aren't influenced by the results. PhysioNet/MIMIC-IV credentialing (CITI training + data use agreement) is in progress as of this date. Once access is granted, Phases 2 onward (cohort extraction, analysis) begin exactly as specified here - no post-hoc changes to eligibility or analysis without an explicitly logged amendment.

---

## 1. Background and rationale

Sepsis resuscitation guidelines (Surviving Sepsis Campaign) mandate an initial fixed-volume fluid challenge before vasopressors are considered, based on a definition of septic shock centered on overt hypotension (Mean Arterial Pressure, MAP < 65 mmHg). This leaves a specific, under-studied patient group poorly addressed by guidelines: patients with **severe tissue hypoperfusion (lactate > 4.0 mmol/L) despite preserved normal blood pressure (MAP > 65 mmHg)** - termed **"Cryptic Shock."** Because these patients look hemodynamically stable, they may not receive timely vasopressor support, even though their lactate indicates significant physiological compromise.

A target trial emulation is the appropriate framework here because a real randomized trial answering this question does not exist, and naive observational comparisons of "early vs. late" treatment are highly susceptible to two specific biases:
- **Confounding by indication** - sicker patients are more likely to receive vasopressors sooner, which can make early treatment look artificially harmful unless properly adjusted for.
- **Immortal time bias** - patients must survive long enough to be classified as "late" treatment, which can make late treatment look artificially better unless time-zero and grace periods are handled correctly.

### 1.1 Directly relevant prior work

This protocol explicitly builds on a **very recent published study** using the same data source and framework:

> Shuaibu II, Radouani S, Hussain Y. "Early Vasopressor Initiation versus Standard Fluid Resuscitation in Normotensive Cryptic Shock: A Target Trial Emulation." *Research Square preprint*, posted February 3, 2026. DOI: 10.21203/rs.3.rs-8747647/v1.

Their study used MIMIC-IV v3.1 (2008–2019), identified 3,558 Cryptic Shock patients (1,245 early-vasopressor, 2,313 standard-care) using Sepsis-3 criteria plus their Cryptic Shock definition, applied a 3-hour grace period at time-zero, adjusted for confounding using IPTW (age, SOFA, lactate as the propensity score model), and reported a large, statistically robust survival benefit for early vasopressor initiation (adjusted HR 0.48, 95% CI 0.40–0.58, p < 0.001), with the benefit persisting across two sensitivity analyses (norepinephrine-only subgroup, and excluding deaths within 6 hours to rule out reverse causation/immortal time artifacts).

A second, related paper - Timing of Core Sepsis Bundle Elements Initiation in Critically Ill Patients: A Multicenter Target Trial Emulation Study (*Clinical Epidemiology*, Dovepress) - independently applies the same TTE-with-IPTW framework to sepsis bundle timing (antibiotics, fluids, vasopressors) across a multicenter cohort, confirming that this exact methodological approach is current, active, and journal-published practice in this field, not a niche technique.

**Why replicate rather than invent a new question from scratch:** Replicating a very recently published, rigorously designed study is a recognized, credible way to build methodological competence - it lets independent verification of the pipeline against a known, published result before extending into original territory. It also means every design decision in Sections 2–5 below is defensible by citation, not by guesswork, which matters given I do not yet have clinical training to independently justify novel threshold choices.

---

## 2. Objective

**Primary objective (replication):** Independently reproduce the Shuaibu et al. (2026) target trial emulation pipeline - cohort definition, IPTW-based confounding adjustment, and Cox survival analysis - using MIMIC-IV, and compare my point estimate and covariate balance diagnostics against their published Table 1 and Table 2 results as a validation check on my own SQL/R pipeline.

**Secondary objective (extension - original contribution):** Having validated the pipeline, extend the analysis in one of the following directions (final choice to be logged as a protocol amendment once Phase 2 cohort numbers are known):

- **(a) Negative control outcome analysis** - test the same exposure (early vasopressor) against an outcome it should *not* plausibly affect (e.g., an unrelated organ-specific outcome not physiologically linked to shock timing). If the "effect" on the negative control is null, this strengthens confidence that the primary result isn't purely an artifact of residual confounding. This is a recognized gap the original paper does not address, and is explicitly flagged in their own Limitations section ("cannot rule out residual unmeasured confounding").
- **(b) Alternative confounder set** - re-run the IPTW model with an expanded covariate set (e.g., adding admission source, vasopressor dose trajectory, or fluid volume more granularly) to test robustness of the HR to model specification.
- **(c) Subgroup analysis** - test whether the treatment effect is consistent across a clinically meaningful subgroup (e.g., patients with vs. without chronic hypertension, since the original discussion specifically hypothesizes the benefit may relate to correcting "relative hypotension" in chronically hypertensive patients - their own mechanistic hypothesis, untested in the paper).

Option (a) is the current lead choice, since it directly answers a limitation the original authors named themselves, is methodologically well-defined, and produces the clearest "I extended published work with an original analysis" story for interviews.

---

## 3. The seven target trial components

| Component | Specification | Source |
|---|---|---|
| **1. Eligibility criteria** | Adults (>18y) meeting Sepsis-3 criteria, with "Cryptic Shock" = simultaneous (i) normotension, MAP > 65 mmHg, and (ii) tissue hypoperfusion, arterial/venous lactate > 4.0 mmol/L. Exclusions: overt hypotension (MAP < 65) prior to T0, active hemorrhage, DNR order within 24h of admission. | Adapted from Shuaibu et al. 2026 |
| **2. Time zero (T0)** | Timestamp of the first lactate measurement > 4.0 mmol/L that is concurrent with a normal MAP measurement. | Adapted from Shuaibu et al. 2026 |
| **3. Treatment strategies** | Classified by intervention received within a 3-hour grace period following T0: **Early Vasopressors** (continuous infusion of norepinephrine, vasopressin, epinephrine, or phenylephrine initiated) vs. **Standard Care** (IV crystalloids only, no vasopressor in the window). | Adapted from Shuaibu et al. 2026 |
| **4. Assignment procedure** | Not randomized (observational); assignment inferred from EHR-recorded interventions within the grace period, using cloning/censoring/weighting logic to avoid immortal time bias for patients whose classification isn't yet determined. | Standard TTE methodology (Hernán & Robins 2016) |
| **5. Follow-up period** | From T0 to 28 days or death/discharge, whichever comes first. | Adapted from Shuaibu et al. 2026 |
| **6. Outcome** | Primary: 28-day all-cause mortality. Secondary (my addition): ICU length of stay. Extension-dependent: negative control outcome (to be specified once cohort is built - candidate: new-onset non-shock-related dermatologic or ophthalmologic diagnosis code, chosen for biological implausibility of a causal link to vasopressor timing). | Primary adapted from paper; secondary/extension original |
| **7. Causal contrast & analysis plan** | Adjusted hazard ratio for 28-day mortality, Early Vasopressors vs. Standard Care, via IPTW-weighted Cox proportional hazards model with robust variance estimation. See Section 4. | Adapted from Shuaibu et al. 2026 |

---

## 4. Statistical analysis plan

1. **Propensity score model:** Multivariable logistic regression predicting treatment (early vasopressor vs. standard care) from baseline confounders measured at/near T0: age, sex, SOFA score, Charlson Comorbidity Index, initial lactate, MAP at T0, and fluid volume administered prior to T0.
2. **Weighting:** Stabilized Inverse Probability of Treatment Weights (IPTW) - weight = 1/PS for treated, 1/(1−PS) for controls, stabilized to the marginal probability of treatment.
3. **Balance diagnostics:** Standardized Mean Differences (SMD) for all covariates before and after weighting; target SMD < 0.1 post-weighting (matches the published paper's threshold, allowing direct comparison of my balance table against theirs).
4. **Outcome analysis:** Weighted Kaplan-Meier survival curves; IPTW-weighted Cox proportional hazards model for the adjusted HR (28-day mortality).
5. **Sensitivity analyses (replicating the paper's):**
   - Restrict treatment definition to norepinephrine only (most common first-line agent).
   - Exclude patients who died within 6 hours of T0 (rules out reverse-causation/immortal-time artifacts).
6. **Extension analysis:** Per the chosen secondary objective in Section 2 - protocol amendment to be logged with date once finalized.
7. **Pre-specified reporting:** Table 1 (baseline characteristics, before/after weighting), Table 2 (unadjusted, multivariable-adjusted, and IPTW-weighted HRs), balance/Love plot, weighted KM curve.

---

## 5. Data source and access

- **Source:** MIMIC-IV v3.1, PhysioNet (Beth Israel Deaconess Medical Center ICU admissions, 2008–2019).
- **Access route:** PhysioNet credentialed access - requires CITI "Data or Specimens Only Research" training completion, credentialing review, and signed Data Use Agreement. **Status: in progress as of 2026-07-13.**
- **Query environment:** Google BigQuery (PhysioNet-hosted MIMIC-IV dataset) for SQL cohort extraction; R for downstream analysis.
- **Relevant tables (planned):** `icustays`, `admissions`, `patients`, `labevents` (lactate), `chartevents` (MAP/vitals), `inputevents` (vasopressor and crystalloid administration), `diagnoses_icd` (Sepsis-3 supporting diagnoses, comorbidities for CCI).

---

## 6. Anticipated limitations

- Single-center data (Beth Israel Deaconess) limits generalizability, same limitation the original paper carries.
- Residual unmeasured confounding is possible despite IPTW - this is precisely why the negative control extension (Option a, Section 2) is a priority if chosen.
- My cohort size and exact patient counts will very likely differ from the published 3,558, due to possible MIMIC-IV version/query differences - this is expected and will be reported transparently, not adjusted to match.
- As a fresher without formal clinical training, threshold choices (MAP > 65, lactate > 4.0) are taken directly from the cited literature rather than independently derived, and this is stated explicitly rather than presented as original clinical judgment.

---

## 7. Literature cited

1. Shuaibu II, Radouani S, Hussain Y. Early Vasopressor Initiation versus Standard Fluid Resuscitation in Normotensive Cryptic Shock: A Target Trial Emulation. *Research Square* [Preprint]. 2026 Feb 3. DOI: 10.21203/rs.3.rs-8747647/v1.
2. [Timing of Core Sepsis Bundle Elements Initiation in Critically Ill Patients: A Multicenter Target Trial Emulation Study](https://www.dovepress.com/timing-of-core-sepsis-bundle-elements-initiation-in-critically-ill-pat-peer-reviewed-fulltext-article-CLEP). *Clinical Epidemiology*, Dovepress.
3. Hernán MA, Robins JM. Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available. *Am J Epidemiol.* 2016;183(8):758-764.
4. Singer M, Deutschman CS, Seymour CW, et al. The Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3). *JAMA.* 2016;315(8):801-810.
5. Johnson AEW, Bulgarelli L, Pollard TJ, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data.* 2023;10(1):1.
6. Austin PC, Stuart EA. Moving towards best practice when using inverse probability of treatment weighting (IPTW) using the propensity score to estimate causal treatment effects in observational studies. *Stat Med.* 2015;34(28):3661-3679.

---

## 8. Amendment log

| Date | Change | Reason |
|---|---|---|
| 2026-07-13 | Protocol v1.0 created | Initial protocol, pre-data-access |
