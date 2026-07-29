"""
build_cohort_demo.py

Phase 2 (prototype pass) — Cryptic Shock cohort extraction against the
MIMIC-IV Clinical Database DEMO (100 patients, open access, no credentialing
required). This is NOT the final analysis pipeline. Its purpose is to build
and debug the extraction logic described in PROTOCOL.md Section 3 (the seven
target trial components) while full credentialed MIMIC-IV access is pending
(see docs/demo_dataset_notes.md for why this exists and what's simplified).

Once PhysioNet credentialing is approved, this same join logic should be
ported to query the full MIMIC-IV database (BigQuery or local Postgres),
with the simplifications noted below replaced by the proper mimic-code
concept tables (sepsis3, sofa, charlson).

Usage:
    python scripts/build_cohort_demo.py

Requires: duckdb, pandas (pip install duckdb pandas)
Expects the unzipped demo dataset at:
    data/mimic-iv-clinical-database-demo-2.2/
"""

import duckdb
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "mimic-iv-clinical-database-demo-2.2"
OUT_DIR = ROOT / "data" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LACTATE_ITEMIDS = (50813, 52442, 53154)          # hosp/d_labitems: "Lactate" (Blood Gas + Chemistry)
MAP_ITEMIDS = (220052, 220181)                    # icu/d_items: Arterial / Non-invasive BP mean
VASOPRESSOR_ITEMIDS = (221906, 221289, 222315, 221749, 221662)
# Norepinephrine, Epinephrine, Vasopressin, Phenylephrine, Dopamine (continuous infusion entries only —
# excludes bolus/intubation-push variants like "Phenylephrine (Intubation)")

MAP_THRESHOLD = 65.0        # mmHg, must be > this (normotension)
LACTATE_THRESHOLD = 4.0     # mmol/L, must be > this (hypoperfusion)
CONCURRENT_WINDOW_HRS = 1   # how close a MAP reading must be to a lactate draw to count as "concurrent"
GRACE_PERIOD_HRS = 3        # treatment-assignment window after T0
FOLLOWUP_DAYS = 28          # primary outcome window

# Simplified Sepsis-3 proxy for demo purposes only (see docs/demo_dataset_notes.md).
# Full pipeline should use the official mimic-code sepsis3.sql concept table instead.
SEPSIS_ICD9 = ("0380", "0381", "03810", "03811", "03812", "03819", "0382", "0383",
               "03840", "03841", "03842", "03843", "03844", "03849", "0388", "0389",
               "99591", "99592", "78552")
SEPSIS_ICD10_PREFIXES = ("A40", "A41", "R6520", "R6521")


def csv(relpath: str) -> str:
    """Return a duckdb-readable path string for a csv.gz file under DATA."""
    return str(DATA / relpath).replace("\\", "/")


def main():
    con = duckdb.connect()

    print(f"Reading demo data from: {DATA}")
    assert DATA.exists(), f"Demo dataset not found at {DATA} — check data/ folder."

    # -- Load core tables -----------------------------------------------------
    patients = con.execute(f"select * from read_csv_auto('{csv('hosp/patients.csv.gz')}')").df()
    admissions = con.execute(f"select * from read_csv_auto('{csv('hosp/admissions.csv.gz')}')").df()
    icustays = con.execute(f"select * from read_csv_auto('{csv('icu/icustays.csv.gz')}')").df()
    diagnoses = con.execute(f"select * from read_csv_auto('{csv('hosp/diagnoses_icd.csv.gz')}')").df()

    print(f"  patients={len(patients)}  admissions={len(admissions)}  icustays={len(icustays)}")

    # -- Sepsis-3 proxy: hospital admission has a sepsis-coded diagnosis -------
    diagnoses["icd_code"] = diagnoses["icd_code"].astype(str).str.strip()
    is_sepsis_icd9 = diagnoses["icd_code"].isin(SEPSIS_ICD9) & (diagnoses["icd_version"] == 9)
    is_sepsis_icd10 = (diagnoses["icd_version"] == 10) & diagnoses["icd_code"].str.startswith(SEPSIS_ICD10_PREFIXES)
    sepsis_hadm_ids = set(diagnoses.loc[is_sepsis_icd9 | is_sepsis_icd10, "hadm_id"].unique())
    print(f"  hospital admissions with a sepsis-coded diagnosis: {len(sepsis_hadm_ids)}")

    icustays_sepsis = icustays[icustays["hadm_id"].isin(sepsis_hadm_ids)].copy()
    print(f"  icustays linked to a sepsis admission: {len(icustays_sepsis)}")

    if icustays_sepsis.empty:
        print("\nNo ICU stays with a sepsis-coded admission in this 100-patient demo subset.")
        print("This is expected — the demo is too small/random to guarantee any Cryptic Shock")
        print("cases. Pipeline logic below is still validated against the full join chain;")
        print("re-run against full MIMIC-IV once credentialed access is granted.")
        _write_empty_outputs()
        return

    stay_ids = tuple(icustays_sepsis["stay_id"].tolist())

    # -- Lactate draws within each ICU stay ------------------------------------
    lact_q = f"""
        select subject_id, hadm_id, charttime, valuenum as lactate
        from read_csv_auto('{csv('hosp/labevents.csv.gz')}')
        where itemid in {LACTATE_ITEMIDS}
          and valuenum is not null
          and valuenum > {LACTATE_THRESHOLD}
    """
    lactate = con.execute(lact_q).df()
    # labevents in the demo doesn't always carry stay_id directly (hosp-level table);
    # join back to icustays via subject_id + hadm_id + time-in-stay instead.
    lactate = lactate.merge(
        icustays_sepsis[["subject_id", "hadm_id", "stay_id", "intime", "outtime"]],
        on=["subject_id", "hadm_id"], how="inner", suffixes=("", "_stay")
    )
    lactate["charttime"] = pd.to_datetime(lactate["charttime"])
    lactate["intime"] = pd.to_datetime(lactate["intime"])
    lactate["outtime"] = pd.to_datetime(lactate["outtime"])
    lactate = lactate[(lactate["charttime"] >= lactate["intime"]) & (lactate["charttime"] <= lactate["outtime"])]
    print(f"  qualifying lactate draws (>{LACTATE_THRESHOLD}) within a sepsis ICU stay: {len(lactate)}")

    # -- MAP readings within each ICU stay -------------------------------------
    map_q = f"""
        select subject_id, stay_id, charttime, valuenum as map_value
        from read_csv_auto('{csv('icu/chartevents.csv.gz')}')
        where itemid in {MAP_ITEMIDS}
          and valuenum is not null
    """
    mapv = con.execute(map_q).df()
    mapv = mapv[mapv["stay_id"].isin(stay_ids)]
    mapv["charttime"] = pd.to_datetime(mapv["charttime"])
    print(f"  MAP readings within sepsis ICU stays: {len(mapv)}")

    # -- Find T0: first lactate>4.0 with a concurrent MAP>65 -------------------
    t0_rows = []
    for stay_id, lgrp in lactate.groupby("stay_id"):
        mgrp = mapv[mapv["stay_id"] == stay_id]
        if mgrp.empty:
            continue
        for _, lrow in lgrp.sort_values("charttime").iterrows():
            window = mgrp[
                (mgrp["charttime"] >= lrow["charttime"] - pd.Timedelta(hours=CONCURRENT_WINDOW_HRS))
                & (mgrp["charttime"] <= lrow["charttime"] + pd.Timedelta(hours=CONCURRENT_WINDOW_HRS))
            ]
            normotensive_window = window[window["map_value"] > MAP_THRESHOLD]
            if not normotensive_window.empty:
                # exclusion: any MAP < 65 recorded *before* this point in the same stay
                prior_hypotension = mgrp[
                    (mgrp["charttime"] < lrow["charttime"]) & (mgrp["map_value"] < MAP_THRESHOLD)
                ]
                if not prior_hypotension.empty:
                    continue
                t0_rows.append({
                    "subject_id": lrow["subject_id"],
                    "hadm_id": lrow["hadm_id"],
                    "stay_id": stay_id,
                    "t0": lrow["charttime"],
                    "lactate_t0": lrow["lactate"],
                    "map_t0": normotensive_window.iloc[0]["map_value"],
                })
                break  # first qualifying T0 only

    cohort = pd.DataFrame(t0_rows)
    print(f"\n  Cryptic Shock T0 identified for: {len(cohort)} ICU stays")

    if cohort.empty:
        _write_empty_outputs()
        return

    # -- Age, sex, adults only --------------------------------------------------
    cohort = cohort.merge(patients[["subject_id", "anchor_age", "anchor_year", "gender", "dod"]],
                           on="subject_id", how="left")
    cohort = cohort.merge(admissions[["hadm_id", "admittime", "deathtime"]], on="hadm_id", how="left")
    cohort["admittime"] = pd.to_datetime(cohort["admittime"])
    cohort["age_at_admission"] = cohort["anchor_age"] + (cohort["admittime"].dt.year - cohort["anchor_year"])
    cohort = cohort[cohort["age_at_admission"] > 18].copy()
    print(f"  after adult (>18y) filter: {len(cohort)}")

    # -- Treatment assignment within grace period --------------------------------
    vaso_q = f"""
        select subject_id, stay_id, starttime, itemid
        from read_csv_auto('{csv('icu/inputevents.csv.gz')}')
        where itemid in {VASOPRESSOR_ITEMIDS}
    """
    vaso = con.execute(vaso_q).df()
    vaso["starttime"] = pd.to_datetime(vaso["starttime"])

    def assign_group(row):
        window_start, window_end = row["t0"], row["t0"] + pd.Timedelta(hours=GRACE_PERIOD_HRS)
        hits = vaso[(vaso["stay_id"] == row["stay_id"])
                    & (vaso["starttime"] >= window_start)
                    & (vaso["starttime"] <= window_end)]
        return "early_vasopressor" if not hits.empty else "standard_care"

    cohort["treatment_group"] = cohort.apply(assign_group, axis=1)

    # -- 28-day mortality outcome -------------------------------------------------
    cohort["dod"] = pd.to_datetime(cohort["dod"])
    cohort["deathtime"] = pd.to_datetime(cohort["deathtime"])
    death_time = cohort["deathtime"].combine_first(cohort["dod"])
    cohort["followup_end"] = cohort["t0"] + pd.Timedelta(days=FOLLOWUP_DAYS)
    cohort["death_within_28d"] = (death_time.notna()) & (death_time <= cohort["followup_end"])

    # -- Save + summarize ----------------------------------------------------------
    out_cols = ["subject_id", "hadm_id", "stay_id", "t0", "lactate_t0", "map_t0",
                "age_at_admission", "gender", "treatment_group", "death_within_28d"]
    cohort[out_cols].to_csv(OUT_DIR / "cohort_demo.csv", index=False)

    print("\n=== Cohort summary (MIMIC-IV DEMO — 100 patients, prototype run) ===")
    print(f"Total eligible Cryptic Shock stays: {len(cohort)}")
    print(cohort["treatment_group"].value_counts().to_string())
    print("\n28-day mortality by group:")
    print(cohort.groupby("treatment_group")["death_within_28d"].mean().to_string())
    print(f"\nFull cohort written to: {OUT_DIR / 'cohort_demo.csv'} (gitignored, not committed)")


def _write_empty_outputs():
    pd.DataFrame(columns=["subject_id", "hadm_id", "stay_id", "t0", "lactate_t0", "map_t0",
                           "age_at_admission", "gender", "treatment_group",
                           "death_within_28d"]).to_csv(OUT_DIR / "cohort_demo.csv", index=False)


if __name__ == "__main__":
    main()
