#---------------------------------------- ========================
######### OASIS-3 FreeSurfer + UDS Data Preparation
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(janitor)

source("R/utilities.R")

ensure_dir("oasis_application/outputs")

# DATA OBJECTS EXPECTED IN THE R ENVIRONMENT
# OASIS3_Freesurfer_output
# OASIS3_UDSb4_cdr
# OASIS3_UDSa1_participant_demo

if (!exists("OASIS3_Freesurfer_output")) {
  stop("OASIS3_Freesurfer_output is not loaded.")
}

if (!exists("OASIS3_UDSb4_cdr")) {
  stop("OASIS3_UDSb4_cdr is not loaded.")
}

if (!exists("OASIS3_UDSa1_participant_demo")) {
  stop("OASIS3_UDSa1_participant_demo is not loaded.")
}

neuro_data  <- OASIS3_Freesurfer_output
neuro_udsb4 <- OASIS3_UDSb4_cdr
neuro_demo  <- OASIS3_UDSa1_participant_demo

### CLEAN COLUMN NAMES

fs_raw <- neuro_data %>%
  clean_names()

b4_raw <- neuro_udsb4 %>%
  clean_names()

demo_raw <- neuro_demo %>%
  clean_names()

cat("FreeSurfer columns:\n")
print(names(fs_raw))

cat("UDS B4 columns:\n")
print(names(b4_raw))

cat("Demo columns:\n")
print(names(demo_raw))

### PREPARE FREESURFER DATA

fs <- fs_raw %>%
  mutate(
    subject = subject,
    mr_session = mr_session,
    fs_id = fs_fsdata_id,
    days_from_mri = as.numeric(str_extract(mr_session, "(?<=_d)[0-9]+")),
    time_years = days_from_mri / 365.25,
    fs_version = factor(version),
    fs_qc_status = fs_qc_status
  ) %>%
  filter(
    !is.na(subject),
    !is.na(mr_session),
    !is.na(intra_cranial_vol)
  )

### PREPARE UDS B4 CDR / COGNITIVE DATA

b4 <- b4_raw %>%
  transmute(
    subject = oasisid,
    uds_session = oasis_session_label,
    days_to_visit_b4 = as.numeric(days_to_visit),
    age_b4 = as.numeric(age_at_visit),
    mmse = as.numeric(mmse),
    memory = as.numeric(memory),
    orient = as.numeric(orient),
    judgment = as.numeric(judgment),
    commun = as.numeric(commun),
    homehobb = as.numeric(homehobb),
    perscare = as.numeric(perscare),
    cdrsum = as.numeric(cdrsum),
    cdrtot = as.numeric(cdrtot),
    dx1_code = dx1_code,
    dx1 = dx1
  ) %>%
  filter(!is.na(subject), !is.na(days_to_visit_b4))

### PREPARE DEMOGRAPHIC FILE

demo <- demo_raw %>%
  transmute(
    subject = oasisid,
    demo_session = oasis_session_label,
    days_to_visit_demo = as.numeric(days_to_visit),
    age_demo = as.numeric(age_at_visit),
    residenc = residenc,
    maristat = maristat,
    livsit = livsit,
    independ = independ
  ) %>%
  filter(!is.na(subject), !is.na(days_to_visit_demo))

### MATCH UDS B4 TO MRI SESSION BY CLOSEST DAYS FROM ENTRY

match_b4_to_mri <- fs %>%
  select(subject, mr_session, days_from_mri) %>%
  left_join(b4, by = "subject") %>%
  mutate(
    abs_day_diff_b4 = abs(days_to_visit_b4 - days_from_mri)
  ) %>%
  group_by(subject, mr_session) %>%
  arrange(abs_day_diff_b4, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(abs_day_diff_b4 <= 365) %>%
  select(
    subject,
    mr_session,
    age_b4,
    mmse,
    memory,
    orient,
    judgment,
    commun,
    homehobb,
    perscare,
    cdrsum,
    cdrtot,
    dx1_code,
    dx1,
    days_to_visit_b4,
    abs_day_diff_b4
  )

### MATCH DEMOGRAPHICS TO MRI SESSION BY CLOSEST DAYS FROM ENTRY

match_demo_to_mri <- fs %>%
  select(subject, mr_session, days_from_mri) %>%
  left_join(demo, by = "subject") %>%
  mutate(
    abs_day_diff_demo = abs(days_to_visit_demo - days_from_mri)
  ) %>%
  group_by(subject, mr_session) %>%
  arrange(abs_day_diff_demo, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(abs_day_diff_demo <= 365) %>%
  select(
    subject,
    mr_session,
    age_demo,
    residenc,
    maristat,
    livsit,
    independ,
    days_to_visit_demo,
    abs_day_diff_demo
  )

### MERGE FREESURFER + UDS B4 + DEMOGRAPHICS

dat0 <- fs %>%
  left_join(match_b4_to_mri, by = c("subject", "mr_session")) %>%
  left_join(match_demo_to_mri, by = c("subject", "mr_session")) %>%
  mutate(
    age = coalesce(age_b4, age_demo),
    age = as.numeric(age)
  ) %>%
  filter(!is.na(age))

cat("Merged data rows:", nrow(dat0), "\n")
cat("Number of subjects:", n_distinct(dat0$subject), "\n")

### SELECT ALZHEIMER-RELEVANT STRUCTURAL REGIONS

region_vars <- c(
  "left_hippocampus_volume",
  "right_hippocampus_volume",
  "total_hippocampus_volume",
  "left_amygdala_volume",
  "right_amygdala_volume",
  "lh_entorhinal_thickness",
  "rh_entorhinal_thickness",
  "lh_parahippocampal_thickness",
  "rh_parahippocampal_thickness",
  "lh_inferiortemporal_thickness",
  "rh_inferiortemporal_thickness",
  "lh_middletemporal_thickness",
  "rh_middletemporal_thickness",
  "lh_precuneus_thickness",
  "rh_precuneus_thickness",
  "lh_posteriorcingulate_thickness",
  "rh_posteriorcingulate_thickness",
  "lh_temporalpole_thickness",
  "rh_temporalpole_thickness"
)

region_vars <- region_vars[region_vars %in% names(dat0)]

cat("Selected regions:\n")
print(region_vars)

### NORMALIZE VOLUME VARIABLES BY INTRACRANIAL VOLUME

volume_vars <- region_vars[str_detect(region_vars, "_volume$")]
thickness_vars <- region_vars[str_detect(region_vars, "_thickness$")]

dat1 <- dat0 %>%
  mutate(
    across(
      all_of(volume_vars),
      ~ .x / intra_cranial_vol,
      .names = "{.col}_icv"
    )
  )

analysis_regions <- c(
  paste0(volume_vars, "_icv"),
  thickness_vars
)

### CREATE LONG FORMAT DATASET

dat_long <- dat1 %>%
  select(
    subject,
    mr_session,
    days_from_mri,
    time_years,
    age,
    sex,
    fs_version,
    fs_qc_status,
    intra_cranial_vol,
    mmse,
    cdrsum,
    cdrtot,
    dx1,
    residenc,
    maristat,
    livsit,
    independ,
    abs_day_diff_b4,
    abs_day_diff_demo,
    all_of(analysis_regions)
  ) %>%
  pivot_longer(
    cols = all_of(analysis_regions),
    names_to = "region",
    values_to = "value"
  ) %>%
  mutate(
    region = str_replace(region, "_icv$", ""),
    region = factor(region),
    subject = factor(subject),
    fs_version = factor(fs_version),
    sex = factor(sex),
    age_c = age - mean(age, na.rm = TRUE),
    age_z = as.numeric(scale(age)),
    time_z = as.numeric(scale(time_years))
  ) %>%
  filter(!is.na(value))

write_csv(dat_long, "oasis_application/outputs/oasis3_model_ready_long.csv")
saveRDS(dat_long, "oasis_application/outputs/oasis3_model_ready_long.rds")

### CREATE WIDE FORMAT DATASET

dat_wide <- dat_long %>%
  select(
    subject,
    mr_session,
    days_from_mri,
    time_years,
    age,
    sex,
    fs_version,
    mmse,
    cdrsum,
    cdrtot,
    dx1,
    region,
    value
  ) %>%
  pivot_wider(
    names_from = region,
    values_from = value
  )

write_csv(dat_wide, "oasis_application/outputs/oasis3_model_ready_wide.csv")
saveRDS(dat_wide, "oasis_application/outputs/oasis3_model_ready_wide.rds")

### SUBJECT-LEVEL LONGITUDINAL SUMMARY

subject_summary <- dat_wide %>%
  group_by(subject) %>%
  summarise(
    n_visits = n_distinct(mr_session),
    baseline_age = min(age, na.rm = TRUE),
    last_age = max(age, na.rm = TRUE),
    followup_years = max(time_years, na.rm = TRUE) - min(time_years, na.rm = TRUE),
    baseline_mmse = mmse[which.min(age)],
    baseline_cdr = cdrtot[which.min(age)],
    .groups = "drop"
  )

summary_overall <- subject_summary %>%
  summarise(
    n_subjects = n(),
    median_visits = median(n_visits, na.rm = TRUE),
    min_visits = min(n_visits, na.rm = TRUE),
    max_visits = max(n_visits, na.rm = TRUE),
    mean_baseline_age = mean(baseline_age, na.rm = TRUE),
    sd_baseline_age = sd(baseline_age, na.rm = TRUE),
    mean_followup_years = mean(followup_years, na.rm = TRUE),
    sd_followup_years = sd(followup_years, na.rm = TRUE),
    mean_baseline_mmse = mean(baseline_mmse, na.rm = TRUE),
    sd_baseline_mmse = sd(baseline_mmse, na.rm = TRUE)
  )

print(summary_overall)

write_csv(subject_summary, "oasis_application/outputs/oasis3_subject_followup_summary.csv")
write_csv(summary_overall, "oasis_application/outputs/oasis3_overall_summary.csv")

### REGION-LEVEL SUMMARY

region_summary <- dat_long %>%
  group_by(region) %>%
  summarise(
    n_obs = n(),
    n_subjects = n_distinct(subject),
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(region_summary, "oasis_application/outputs/oasis3_region_summary.csv")

### KEEP SUBJECTS WITH AT LEAST TWO MRI VISITS

subjects_two_plus <- dat_wide %>%
  count(subject, name = "n_visits") %>%
  filter(n_visits >= 2)

dat_long_2plus <- dat_long %>%
  semi_join(subjects_two_plus, by = "subject")

dat_wide_2plus <- dat_wide %>%
  semi_join(subjects_two_plus, by = "subject")

write_csv(dat_long_2plus, "oasis_application/outputs/oasis3_model_ready_long_2plus.csv")
saveRDS(dat_long_2plus, "oasis_application/outputs/oasis3_model_ready_long_2plus.rds")

write_csv(dat_wide_2plus, "oasis_application/outputs/oasis3_model_ready_wide_2plus.csv")
saveRDS(dat_wide_2plus, "oasis_application/outputs/oasis3_model_ready_wide_2plus.rds")

#"Subjects with >=2 visits
n_distinct(dat_long_2plus$subject)
#"Rows in long 2+ dataset:
 nrow(dat_long_2plus)

