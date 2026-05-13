#---------------------------------------- ========================
######### OASIS-3 Manuscript and Supplementary Tables
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(readr)

source("R/utilities.R")

ensure_dir("manuscript/tables")
ensure_dir("manuscript/supplementary/tables")
ensure_dir("oasis_application/outputs/tables")

analysis_objects <- readRDS(
  "oasis_application/outputs/model_objects/oasis_analysis_objects.rds"
)

compare_tbl <- analysis_objects$compare_tbl
region_dev_C <- analysis_objects$region_dev_C
subject_dev_C <- analysis_objects$subject_dev_C
fit_C_summary <- analysis_objects$fit_C_summary
key_summary <- analysis_objects$key_summary

### MAIN TABLES

table_model_comparison <- compare_tbl %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

table_top_regions <- region_dev_C %>%
  mutate(region_clean = clean_region_name(region)) %>%
  arrange(desc(tail_prob)) %>%
  slice_head(n = 10) %>%
  select(region_clean, n_obs, mean_z, sd_z, tail_prob) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

table_posterior_key <- key_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write_csv(table_model_comparison, "manuscript/tables/Table1_Model_Comparison.csv")
write_csv(table_posterior_key, "manuscript/tables/Table2_Key_Posterior_Parameters.csv")
write_csv(table_top_regions, "manuscript/tables/Table3_Top_Regions_Extreme_Deviations.csv")

### SUPPLEMENTARY TABLES

write_csv(
  region_dev_C %>%
    mutate(region_clean = clean_region_name(region)) %>%
    select(region_clean, n_obs, mean_z, sd_z, tail_prob) %>%
    arrange(desc(tail_prob)),
  "manuscript/supplementary/tables/TableS1_All_Region_Deviation_Summaries.csv"
)

write_csv(
  subject_dev_C %>%
    arrange(desc(mean_abs_z)),
  "manuscript/supplementary/tables/TableS2_All_Subject_Deviation_Summaries.csv"
)

write_csv(
  key_summary,
  "manuscript/supplementary/tables/TableS3_Key_Posterior_Parameter_Summary.csv"
)

write_csv(
  fit_C_summary,
  "manuscript/supplementary/tables/TableS4_Full_Stan_Parameter_Summary.csv"
)

