# ---------------------------------------------------------------------------
# pipeline_demo_broadened.R
#
# PIPELINE DEMONSTRATION ONLY - NOT A SCIENTIFIC RESULT.
#
# Same broadened-cohort pipeline as the Python version (all sepsis ICU stays
# in the MIMIC-IV Demo, split by "received a vasopressor" vs. "did not" -
# NOT the strict Cryptic Shock definition). Exists to build/debug the
# propensity score, IPTW, balance diagnostics, and Cox model code ahead of
# full credentialed MIMIC-IV access. Do not report any output as a finding.
#
# Requires: readr, dplyr, lubridate, survival, ggplot2, tidyr, stringr
# ---------------------------------------------------------------------------

library(readr)
library(dplyr)
library(lubridate)
library(survival)
library(ggplot2)
library(tidyr)
library(stringr)

ROOT <- getwd()
# If not running in RStudio, just set this manually instead:
# ROOT <- "path/to/TTE_Cryptic_Shock_Project"

DATA_DIR   <- file.path(ROOT, "data", "mimic-iv-clinical-database-demo-2.2")
OUT_DIR    <- file.path(ROOT, "data", "processed")
REPORT_DIR <- file.path(ROOT, "report")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

LACTATE_ITEMIDS     <- c(50813, 52442, 53154)
MAP_ITEMIDS         <- c(220052, 220181)
VASOPRESSOR_ITEMIDS <- c(221906, 221289, 222315, 221749, 221662)
FOLLOWUP_DAYS       <- 28

SEPSIS_ICD9 <- c("0380","0381","03810","03811","03812","03819","0382","0383",
                 "03840","03841","03842","03843","03844","03849","0388","0389",
                 "99591","99592","78552")
SEPSIS_ICD10_PREFIXES <- c("A40", "A41", "R6520", "R6521")

read_gz <- function(relpath) {
  read_csv(file.path(DATA_DIR, relpath), show_col_types = FALSE)
}

# readr auto-detects date/datetime columns on read, so by the time a column
# reaches our code it may already be POSIXct/Date rather than character.
# Calling ymd_hms() directly on an already-parsed column can silently fail
# to parse (this is what caused "270 failed to parse" earlier). Wrapping in
# as.character() first, then using parse_date_time() with multiple possible
# formats, makes this robust regardless of what readr already did.
parse_dt <- function(x) {
  parse_date_time(as.character(x), orders = c("Ymd HMS", "Ymd HM", "Ymd"), quiet = TRUE)
}

# ---------------------------------------------------------------------------
# 1. Build the broadened cohort
# ---------------------------------------------------------------------------

build_cohort <- function() {
  patients   <- read_gz("hosp/patients.csv.gz")
  admissions <- read_gz("hosp/admissions.csv.gz")
  icustays   <- read_gz("icu/icustays.csv.gz")
  diagnoses  <- read_gz("hosp/diagnoses_icd.csv.gz") %>%
    mutate(icd_code = str_trim(as.character(icd_code)))
  
  is_sepsis9  <- diagnoses$icd_code %in% SEPSIS_ICD9 & diagnoses$icd_version == 9
  is_sepsis10 <- diagnoses$icd_version == 10 &
    str_starts(diagnoses$icd_code, paste(SEPSIS_ICD10_PREFIXES, collapse = "|"))
  sepsis_hadm_ids <- unique(diagnoses$hadm_id[is_sepsis9 | is_sepsis10])
  
  comorbidity_count <- diagnoses %>%
    group_by(hadm_id) %>%
    summarise(comorbidity_count = n_distinct(icd_code), .groups = "drop")
  
  cohort <- icustays %>%
    filter(hadm_id %in% sepsis_hadm_ids) %>%
    left_join(patients %>% select(subject_id, anchor_age, anchor_year, gender, dod),
              by = "subject_id") %>%
    left_join(admissions %>% select(hadm_id, admittime, deathtime), by = "hadm_id") %>%
    left_join(comorbidity_count, by = "hadm_id") %>%
    mutate(
      intime    = parse_dt(intime),
      outtime   = parse_dt(outtime),
      admittime = parse_dt(admittime),
      age_at_admission = anchor_age + (year(admittime) - anchor_year)
    ) %>%
    filter(age_at_admission > 18)
  
  # -- treatment: any vasopressor during the ICU stay -------------------
  vaso <- read_gz("icu/inputevents.csv.gz") %>%
    filter(itemid %in% VASOPRESSOR_ITEMIDS)
  treated_stays <- unique(vaso$stay_id)
  cohort <- cohort %>% mutate(treatment = as.integer(stay_id %in% treated_stays))
  
  # -- baseline covariates: first lactate / first MAP during the stay ---
  lact <- read_gz("hosp/labevents.csv.gz") %>%
    filter(itemid %in% LACTATE_ITEMIDS, !is.na(valuenum)) %>%
    rename(lactate = valuenum) %>%
    mutate(charttime = parse_dt(charttime)) %>%
    inner_join(cohort %>% select(subject_id, hadm_id, stay_id, intime, outtime),
               by = c("subject_id", "hadm_id"), relationship = "many-to-many") %>%
    filter(charttime >= intime, charttime <= outtime) %>%
    arrange(charttime) %>%
    group_by(stay_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(stay_id, initial_lactate = lactate)
  
  mapv <- read_gz("icu/chartevents.csv.gz") %>%
    filter(itemid %in% MAP_ITEMIDS, !is.na(valuenum)) %>%
    rename(map_value = valuenum) %>%
    mutate(charttime = parse_dt(charttime)) %>%
    arrange(charttime) %>%
    group_by(stay_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(stay_id, initial_map = map_value)
  
  cohort <- cohort %>%
    left_join(lact, by = "stay_id") %>%
    left_join(mapv, by = "stay_id")
  
  # -- outcome: 28-day mortality from ICU intime -------------------------
  cohort <- cohort %>%
    mutate(
      dod          = parse_dt(dod),
      deathtime    = parse_dt(deathtime),
      death_time   = coalesce(deathtime, dod),
      followup_end = intime + days(FOLLOWUP_DAYS),
      died_in_window = !is.na(death_time) & death_time <= followup_end,
      event_time   = if_else(died_in_window, death_time, followup_end),
      time_days    = as.numeric(difftime(event_time, intime, units = "days")),
      event_28d    = as.integer(died_in_window),
      gender_male  = as.integer(gender == "M")
    ) %>%
    select(subject_id, hadm_id, stay_id, age_at_admission, gender_male,
           initial_lactate, initial_map, comorbidity_count,
           treatment, time_days, event_28d) %>%
    drop_na(initial_lactate, initial_map)
  
  cohort
}

# ---------------------------------------------------------------------------
# 2. Propensity score model + stabilized IPTW weights
# ---------------------------------------------------------------------------

COVARIATES <- c("age_at_admission", "gender_male", "initial_lactate",
                "initial_map", "comorbidity_count")

fit_propensity_and_weights <- function(cohort) {
  form <- as.formula(paste("treatment ~", paste(COVARIATES, collapse = " + ")))
  ps_model <- glm(form, data = cohort, family = binomial())
  cohort$ps <- predict(ps_model, type = "response")
  
  marginal_p <- mean(cohort$treatment)
  cohort$iptw <- ifelse(
    cohort$treatment == 1,
    marginal_p / cohort$ps,
    (1 - marginal_p) / (1 - cohort$ps)
  )
  list(cohort = cohort, model = ps_model)
}

# ---------------------------------------------------------------------------
# 3. Balance diagnostics: SMD before/after weighting
# ---------------------------------------------------------------------------

weighted_var <- function(x, w, mean_x) {
  sum(w * (x - mean_x)^2) / sum(w)
}

smd_continuous <- function(x, treat, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  t_idx <- treat == 1
  c_idx <- treat == 0
  m1 <- weighted.mean(x[t_idx], w[t_idx])
  m0 <- weighted.mean(x[c_idx], w[c_idx])
  v1 <- weighted_var(x[t_idx], w[t_idx], m1)
  v0 <- weighted_var(x[c_idx], w[c_idx], m0)
  (m1 - m0) / sqrt((v1 + v0) / 2)
}

balance_table <- function(cohort) {
  rows <- lapply(COVARIATES, function(cov) {
    unw <- smd_continuous(cohort[[cov]], cohort$treatment)
    w   <- smd_continuous(cohort[[cov]], cohort$treatment, cohort$iptw)
    data.frame(covariate = cov, smd_unweighted = unw, smd_weighted = w)
  })
  bind_rows(rows)
}

love_plot <- function(balance_df, path) {
  plot_df <- balance_df %>%
    pivot_longer(cols = c(smd_unweighted, smd_weighted),
                 names_to = "type", values_to = "smd") %>%
    mutate(smd = abs(smd),
           type = recode(type, smd_unweighted = "Unweighted", smd_weighted = "IPTW-weighted"))
  
  p <- ggplot(plot_df, aes(x = smd, y = covariate, color = type)) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray40") +
    labs(x = "Absolute Standardized Mean Difference", y = NULL, color = NULL,
         title = "Covariate balance before/after IPTW (pipeline demo)") +
    theme_minimal()
  
  ggsave(path, p, width = 6, height = 4, dpi = 150)
}

# ---------------------------------------------------------------------------
# 4. Weighted Cox model + weighted KM curve
# ---------------------------------------------------------------------------

fit_weighted_cox <- function(cohort) {
  coxph(Surv(time_days, event_28d) ~ treatment, data = cohort,
        weights = iptw, robust = TRUE)
}

weighted_km_plot <- function(cohort, path) {
  fit <- survfit(Surv(time_days, event_28d) ~ treatment, data = cohort, weights = iptw)
  
  png(path, width = 6, height = 4, units = "in", res = 150)
  plot(fit, col = c("red", "blue"), lwd = 2,
       xlab = "Days from ICU admission", ylab = "Survival probability",
       main = "IPTW-weighted Kaplan-Meier (pipeline demo - not a result)")
  legend("bottomleft", legend = c("Standard care", "Early vasopressor"),
         col = c("red", "blue"), lwd = 2)
  dev.off()
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  cohort <- build_cohort()
  cat(sprintf("Broadened demo cohort: n=%d (treated=%d, control=%d)\n",
              nrow(cohort), sum(cohort$treatment), sum(cohort$treatment == 0)))
  
  if (n_distinct(cohort$treatment) < 2 || nrow(cohort) < 10) {
    cat("Not enough patients/variation in this demo pull to fit the pipeline. Stopping.\n")
    return(invisible(NULL))
  }
  
  res <- fit_propensity_and_weights(cohort)
  cohort <- res$cohort
  cat("\nPropensity score model (logistic regression) summary:\n")
  print(summary(res$model))
  
  balance_df <- balance_table(cohort)
  cat("\nBalance table (SMD before/after IPTW):\n")
  print(balance_df)
  love_plot(balance_df, file.path(REPORT_DIR, "love_plot_demo_R.png"))
  
  cox_model <- fit_weighted_cox(cohort)
  cat("\nWeighted Cox model (treatment effect on 28-day mortality):\n")
  print(summary(cox_model))
  
  weighted_km_plot(cohort, file.path(REPORT_DIR, "km_curve_demo_R.png"))
  
  write_csv(cohort, file.path(OUT_DIR, "broadened_cohort_demo_R.csv"))
  write_csv(balance_df, file.path(REPORT_DIR, "table1_balance_demo_R.csv"))
  write_csv(as.data.frame(summary(cox_model)$coefficients),
            file.path(REPORT_DIR, "table2_cox_demo_R.csv"))
  
  cat("\nREMINDER: this is a pipeline demonstration on a loosened, ~20-patient\n")
  cat("subgroup of the 100-patient MIMIC-IV Demo. None of these numbers are a\n")
  cat("scientific finding - they exist only to prove the code runs correctly\n")
  cat("before being pointed at the real, credentialed MIMIC-IV cohort.\n")
}

main()
