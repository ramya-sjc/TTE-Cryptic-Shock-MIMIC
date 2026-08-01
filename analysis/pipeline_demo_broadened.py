"""
pipeline_demo_broadened.py

PIPELINE DEMONSTRATION ONLY — NOT A SCIENTIFIC RESULT.

Runs the full downstream statistical pipeline (propensity score model, IPTW
weighting, balance diagnostics, weighted Cox model, weighted KM curve)
against a BROADENED cohort from the open MIMIC-IV Demo (100 patients):
all sepsis ICU stays, split by "received any vasopressor during the stay"
vs. "did not" — NOT the strict Cryptic Shock definition in PROTOCOL.md.

Why broadened: the strict Cryptic Shock T0 definition (see
scripts/build_cohort_demo.py) found only 1 eligible patient in the 100-patient
demo — too few for any weighting/regression method to run meaningfully. This
script exists purely to build and debug the code for every remaining phase
(Sections 4.1-4.4 of PROTOCOL.md) so it's ready to point at the real,
full-scale MIMIC-IV cohort once PhysioNet credentialing is approved.

Do not cite, report, or present any number this script outputs as a finding.
Sample size (~20-something patients) is far too small for a valid estimate.

Requires: duckdb, pandas, numpy, statsmodels, lifelines, matplotlib
"""

import duckdb
import numpy as np
import pandas as pd
import statsmodels.api as sm
from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "mimic-iv-clinical-database-demo-2.2"
OUT_DIR = ROOT / "data" / "processed"
REPORT_DIR = ROOT / "report"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

LACTATE_ITEMIDS = (50813, 52442, 53154)
MAP_ITEMIDS = (220052, 220181)
VASOPRESSOR_ITEMIDS = (221906, 221289, 222315, 221749, 221662)
FOLLOWUP_DAYS = 28

SEPSIS_ICD9 = ("0380", "0381", "03810", "03811", "03812", "03819", "0382", "0383",
               "03840", "03841", "03842", "03843", "03844", "03849", "0388", "0389",
               "99591", "99592", "78552")
SEPSIS_ICD10_PREFIXES = ("A40", "A41", "R6520", "R6521")


def csv(relpath: str) -> str:
    return str(DATA / relpath).replace("\\", "/")


# ---------------------------------------------------------------------------
# 1. Build the broadened cohort
# ---------------------------------------------------------------------------

def build_cohort(con):
    patients = con.execute(f"select * from read_csv_auto('{csv('hosp/patients.csv.gz')}')").df()
    admissions = con.execute(f"select * from read_csv_auto('{csv('hosp/admissions.csv.gz')}')").df()
    icustays = con.execute(f"select * from read_csv_auto('{csv('icu/icustays.csv.gz')}')").df()
    diagnoses = con.execute(f"select * from read_csv_auto('{csv('hosp/diagnoses_icd.csv.gz')}')").df()

    diagnoses["icd_code"] = diagnoses["icd_code"].astype(str).str.strip()
    is_sepsis9 = diagnoses["icd_code"].isin(SEPSIS_ICD9) & (diagnoses["icd_version"] == 9)
    is_sepsis10 = (diagnoses["icd_version"] == 10) & diagnoses["icd_code"].str.startswith(SEPSIS_ICD10_PREFIXES)
    sepsis_hadm_ids = set(diagnoses.loc[is_sepsis9 | is_sepsis10, "hadm_id"].unique())

    comorbidity_count = diagnoses.groupby("hadm_id")["icd_code"].nunique().rename("comorbidity_count")

    cohort = icustays[icustays["hadm_id"].isin(sepsis_hadm_ids)].copy()
    cohort = cohort.merge(patients[["subject_id", "anchor_age", "anchor_year", "gender", "dod"]],
                           on="subject_id", how="left")
    cohort = cohort.merge(admissions[["hadm_id", "admittime", "deathtime"]], on="hadm_id", how="left")
    cohort = cohort.merge(comorbidity_count, on="hadm_id", how="left")

    cohort["intime"] = pd.to_datetime(cohort["intime"])
    cohort["outtime"] = pd.to_datetime(cohort["outtime"])
    cohort["admittime"] = pd.to_datetime(cohort["admittime"])
    cohort["age_at_admission"] = cohort["anchor_age"] + (cohort["admittime"].dt.year - cohort["anchor_year"])
    cohort = cohort[cohort["age_at_admission"] > 18].copy()

    # -- treatment: any vasopressor during the ICU stay --------------------
    vaso = con.execute(f"""
        select stay_id, starttime from read_csv_auto('{csv('icu/inputevents.csv.gz')}')
        where itemid in {VASOPRESSOR_ITEMIDS}
    """).df()
    treated_stays = set(vaso["stay_id"].unique())
    cohort["treatment"] = cohort["stay_id"].isin(treated_stays).astype(int)

    # -- baseline covariates: first lactate / first MAP during the stay ----
    lact = con.execute(f"""
        select subject_id, hadm_id, charttime, valuenum as lactate
        from read_csv_auto('{csv('hosp/labevents.csv.gz')}')
        where itemid in {LACTATE_ITEMIDS} and valuenum is not null
    """).df()
    lact["charttime"] = pd.to_datetime(lact["charttime"])
    lact = lact.merge(cohort[["subject_id", "hadm_id", "stay_id", "intime", "outtime"]],
                       on=["subject_id", "hadm_id"], how="inner")
    lact = lact[(lact["charttime"] >= lact["intime"]) & (lact["charttime"] <= lact["outtime"])]
    first_lactate = lact.sort_values("charttime").groupby("stay_id").first()["lactate"].rename("initial_lactate")

    mapv = con.execute(f"""
        select stay_id, charttime, valuenum as map_value
        from read_csv_auto('{csv('icu/chartevents.csv.gz')}')
        where itemid in {MAP_ITEMIDS} and valuenum is not null
    """).df()
    mapv["charttime"] = pd.to_datetime(mapv["charttime"])
    first_map = mapv.sort_values("charttime").groupby("stay_id").first()["map_value"].rename("initial_map")

    cohort = cohort.merge(first_lactate, on="stay_id", how="left")
    cohort = cohort.merge(first_map, on="stay_id", how="left")

    # -- outcome: 28-day mortality from ICU intime --------------------------
    cohort["dod"] = pd.to_datetime(cohort["dod"])
    cohort["deathtime"] = pd.to_datetime(cohort["deathtime"])
    death_time = cohort["deathtime"].combine_first(cohort["dod"])
    followup_end = cohort["intime"] + pd.Timedelta(days=FOLLOWUP_DAYS)
    died_in_window = death_time.notna() & (death_time <= followup_end)
    event_time = death_time.where(died_in_window, followup_end)
    cohort["time_days"] = (event_time - cohort["intime"]).dt.total_seconds() / 86400
    cohort["event_28d"] = died_in_window.astype(int)

    keep_cols = ["subject_id", "hadm_id", "stay_id", "age_at_admission", "gender",
                 "initial_lactate", "initial_map", "comorbidity_count",
                 "treatment", "time_days", "event_28d"]
    cohort = cohort[keep_cols].dropna(subset=["initial_lactate", "initial_map"]).reset_index(drop=True)
    cohort["gender_male"] = (cohort["gender"] == "M").astype(int)
    return cohort


# ---------------------------------------------------------------------------
# 2. Propensity score model + stabilized IPTW weights
# ---------------------------------------------------------------------------

COVARIATES = ["age_at_admission", "gender_male", "initial_lactate", "initial_map", "comorbidity_count"]


def fit_propensity_and_weights(cohort):
    X = sm.add_constant(cohort[COVARIATES])
    y = cohort["treatment"]
    model = sm.Logit(y, X).fit(disp=0)
    cohort["ps"] = model.predict(X)

    marginal_p = cohort["treatment"].mean()
    cohort["iptw"] = np.where(
        cohort["treatment"] == 1,
        marginal_p / cohort["ps"],
        (1 - marginal_p) / (1 - cohort["ps"]),
    )
    return cohort, model


# ---------------------------------------------------------------------------
# 3. Balance diagnostics: SMD before/after weighting
# ---------------------------------------------------------------------------

def smd_continuous(x, treat, weights=None):
    if weights is None:
        weights = np.ones(len(x))
    t, c = treat == 1, treat == 0
    m1 = np.average(x[t], weights=weights[t])
    m0 = np.average(x[c], weights=weights[c])
    v1 = np.average((x[t] - m1) ** 2, weights=weights[t])
    v0 = np.average((x[c] - m0) ** 2, weights=weights[c])
    return (m1 - m0) / np.sqrt((v1 + v0) / 2)


def balance_table(cohort):
    rows = []
    for cov in COVARIATES:
        unw = smd_continuous(cohort[cov].values, cohort["treatment"].values)
        w = smd_continuous(cohort[cov].values, cohort["treatment"].values, cohort["iptw"].values)
        rows.append({"covariate": cov, "smd_unweighted": unw, "smd_weighted": w})
    return pd.DataFrame(rows)


def love_plot(balance_df, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    y = np.arange(len(balance_df))
    ax.scatter(balance_df["smd_unweighted"].abs(), y, label="Unweighted", color="tab:red")
    ax.scatter(balance_df["smd_weighted"].abs(), y, label="IPTW-weighted", color="tab:blue")
    ax.axvline(0.1, linestyle="--", color="gray", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(balance_df["covariate"])
    ax.set_xlabel("Absolute Standardized Mean Difference")
    ax.set_title("Covariate balance before/after IPTW (pipeline demo)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 4. Weighted Cox model + weighted KM curve
# ---------------------------------------------------------------------------

def fit_weighted_cox(cohort):
    df = cohort[["time_days", "event_28d", "treatment", "iptw"]].copy()
    cph = CoxPHFitter()
    cph.fit(df, duration_col="time_days", event_col="event_28d",
            weights_col="iptw", robust=True)
    return cph


def weighted_km_plot(cohort, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    kmf = KaplanMeierFitter()
    for label, grp in cohort.groupby("treatment"):
        name = "Early vasopressor" if label == 1 else "Standard care"
        kmf.fit(grp["time_days"], event_observed=grp["event_28d"],
                weights=grp["iptw"], label=name)
        kmf.plot_survival_function(ax=ax)
    ax.set_xlabel("Days from ICU admission")
    ax.set_ylabel("Survival probability")
    ax.set_title("IPTW-weighted Kaplan-Meier (pipeline demo — not a result)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    con = duckdb.connect()
    cohort = build_cohort(con)
    print(f"Broadened demo cohort: n={len(cohort)}  "
          f"(treated={cohort['treatment'].sum()}, control={(cohort['treatment']==0).sum()})")

    if cohort["treatment"].nunique() < 2 or len(cohort) < 10:
        print("Not enough patients/variation in this demo pull to fit the pipeline. Stopping.")
        return

    cohort, ps_model = fit_propensity_and_weights(cohort)
    print("\nPropensity score model (logistic regression) summary:")
    print(ps_model.summary2().tables[1])

    balance_df = balance_table(cohort)
    print("\nBalance table (SMD before/after IPTW):")
    print(balance_df.to_string(index=False))
    love_plot(balance_df, REPORT_DIR / "love_plot_demo.png")

    cph = fit_weighted_cox(cohort)
    print("\nWeighted Cox model (treatment effect on 28-day mortality):")
    print(cph.summary[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]])

    weighted_km_plot(cohort, REPORT_DIR / "km_curve_demo.png")

    cohort.to_csv(OUT_DIR / "broadened_cohort_demo.csv", index=False)
    balance_df.to_csv(REPORT_DIR / "table1_balance_demo.csv", index=False)
    cph.summary.to_csv(REPORT_DIR / "table2_cox_demo.csv")

    print(f"\nOutputs written to {REPORT_DIR} and {OUT_DIR} (gitignored data; plots/tables are fine to commit).")
    print("\nREMINDER: this is a pipeline demonstration on a loosened, ~20-patient")
    print("subgroup of the 100-patient MIMIC-IV Demo. None of these numbers are a")
    print("scientific finding — they exist only to prove the code runs correctly")
    print("before being pointed at the real, credentialed MIMIC-IV cohort.")


if __name__ == "__main__":
    main()
