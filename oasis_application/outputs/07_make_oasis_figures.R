#---------------------------------------- ========================
######### OASIS-3 Manuscript and Supplementary Figures
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)
library(forcats)
library(patchwork)
library(viridis)
library(ggseg)
library(bayesplot)

source("R/utilities.R")

ensure_dir("manuscript/figures")
ensure_dir("manuscript/supplementary/figures")
ensure_dir("manuscript/supplementary/tables")

analysis_objects <- readRDS(
  "oasis_application/outputs/model_objects/oasis_analysis_objects.rds"
)

fit_C_stan <- readRDS(
  "oasis_application/outputs/model_objects/fit_C_subject_specific_spatial_stan.rds"
)

dat_stan <- analysis_objects$dat_stan
compare_tbl <- analysis_objects$compare_tbl
region_dev_C <- analysis_objects$region_dev_C
subject_dev_C <- analysis_objects$subject_dev_C
fit_C_summary <- analysis_objects$fit_C_summary
key_summary <- analysis_objects$key_summary

### MAIN FIGURE 1: MODEL COMPARISON

compare_long <- compare_tbl %>%
  select(model, rmse, mean_abs_deviation) %>%
  pivot_longer(
    cols = c(rmse, mean_abs_deviation),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      rmse = "Root mean squared error",
      mean_abs_deviation = "Mean absolute deviation"
    ),
    model = factor(
      model,
      levels = c(
        "Independent cross-sectional model",
        "Longitudinal non-spatial model",
        "Bayesian subject-specific spatial model"
      )
    )
  )

p_model_compare <- ggplot(compare_long, aes(x = model, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), color = "black", width = 0.65) +
  coord_flip() +
  labs(
    subtitle = "Lower values indicate better subject-level normative fit",
    x = NULL,
    y = "Error",
    fill = NULL
  ) +
  theme_pub(13)

ggsave(
  "manuscript/figures/Figure1_Model_Comparison_RMSE_MAD.pdf",
  p_model_compare,
  width = 8.5,
  height = 5.2,
  bg = "white"
)

ggsave(
  "manuscript/figures/Figure1_Model_Comparison_RMSE_MAD.png",
  p_model_compare,
  width = 8.5,
  height = 5.2,
  dpi = 500,
  bg = "white"
)

### MAIN FIGURE 2: DK ATLAS EXTREME DEVIATION MAP

dk_regions <- ggseg.formats::atlas_regions(dk())

atlas_df <- region_dev_C %>%
  mutate(
    original_region = as.character(region),
    hemi = case_when(
      str_detect(original_region, "^lh_") ~ "left",
      str_detect(original_region, "^rh_") ~ "right",
      TRUE ~ NA_character_
    ),
    region = original_region %>%
      str_replace("^lh_", "") %>%
      str_replace("^rh_", "") %>%
      str_replace("_thickness$", "") %>%
      str_replace("_volume$", "") %>%
      str_replace("inferiortemporal", "inferior temporal") %>%
      str_replace("middletemporal", "middle temporal") %>%
      str_replace("posteriorcingulate", "posterior cingulate") %>%
      str_replace("temporalpole", "temporal pole"),
    region_display = case_when(
      region == "entorhinal" ~ "Entorhinal cortex",
      region == "parahippocampal" ~ "Parahippocampal cortex",
      region == "posterior cingulate" ~ "Posterior cingulate",
      region == "precuneus" ~ "Precuneus",
      region == "inferior temporal" ~ "Inferior temporal",
      region == "middle temporal" ~ "Middle temporal",
      region == "temporal pole" ~ "Temporal pole",
      TRUE ~ str_to_title(region)
    ),
    hemi_display = case_when(
      hemi == "left" ~ "Left",
      hemi == "right" ~ "Right",
      TRUE ~ hemi
    ),
    region_hemi = paste(hemi_display, region_display)
  ) %>%
  filter(!is.na(hemi), region %in% dk_regions)

write_csv(
  atlas_df %>%
    select(original_region, hemi, region, region_hemi, mean_z, sd_z, tail_prob),
  "manuscript/supplementary/tables/TableS_DK_Mapped_Cortical_Regions.csv"
)

tail_limits <- range(atlas_df$tail_prob, na.rm = TRUE)
mean_lim <- max(abs(atlas_df$mean_z), na.rm = TRUE)

p_brain_tail <- ggplot(atlas_df) +
  geom_brain(
    atlas = dk(),
    position = position_brain(hemi ~ view),
    aes(fill = tail_prob)
  ) +
  scale_fill_viridis_c(
    option = "cividis",
    direction = -1,
    limits = tail_limits,
    labels = scales::number_format(accuracy = 0.001),
    na.value = "grey85"
  ) +
  labs(
    subtitle = "Probability of |z| > 1.96 under the Bayesian spatial normative model",
    fill = expression(Pr("|z| > 1.96"))
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

p_region_tail_panel <- atlas_df %>%
  arrange(desc(tail_prob)) %>%
  ggplot(aes(x = fct_reorder(region_hemi, tail_prob), y = tail_prob)) +
  geom_col(color = "black", width = 0.7) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(Pr("|z| > 1.96"))
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

p_brain_tail_labeled <- p_brain_tail + p_region_tail_panel +
  plot_layout(widths = c(2.2, 1))

ggsave(
  "manuscript/figures/Figure2_Cortical_Extreme_Deviations.pdf",
  p_brain_tail_labeled,
  width = 15,
  height = 7,
  bg = "white"
)

ggsave(
  "manuscript/figures/Figure2_Cortical_Extreme_Deviations.png",
  p_brain_tail_labeled,
  width = 15,
  height = 7,
  dpi = 600,
  bg = "white"
)

### MAIN FIGURE 3: LONGITUDINAL TRAJECTORIES

key_regions <- c(
  "left_hippocampus_volume",
  "right_hippocampus_volume",
  "left_amygdala_volume",
  "right_amygdala_volume",
  "lh_entorhinal_thickness",
  "rh_entorhinal_thickness",
  "lh_parahippocampal_thickness",
  "rh_parahippocampal_thickness"
)

key_regions <- key_regions[key_regions %in% levels(dat_stan$region)]

dat_traj <- dat_stan %>%
  filter(region %in% key_regions) %>%
  mutate(region_clean = clean_region_name(region))

p_dev_traj <- ggplot(dat_traj, aes(x = age, y = z_C, group = subject)) +
  geom_line(alpha = 0.10, linewidth = 0.3) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, linewidth = 0.9) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ region_clean, scales = "free_y", ncol = 2) +
  labs(
    subtitle = "Subject-specific deviations from the Bayesian spatial normative model",
    x = "Age at MRI visit",
    y = "Standardized deviation"
  ) +
  theme_pub(13)

ggsave(
  "manuscript/figures/Figure3_Longitudinal_Deviation_Trajectories.pdf",
  p_dev_traj,
  width = 10.5,
  height = 8,
  bg = "white"
)

ggsave(
  "manuscript/figures/Figure3_Longitudinal_Deviation_Trajectories.png",
  p_dev_traj,
  width = 10.5,
  height = 8,
  dpi = 500,
  bg = "white"
)

### MAIN FIGURE 4: SUBJECT CASE STUDIES

subject_rank <- subject_dev_C %>%
  filter(!is.na(mean_abs_z), n_visits >= 2) %>%
  arrange(mean_abs_z)

abnormal_subject <- subject_rank %>%
  arrange(desc(mean_abs_z)) %>%
  slice(1) %>%
  pull(subject)

normal_subject <- subject_rank %>%
  filter(
    n_visits >= 3,
    mean_abs_z > 0.30,
    mean_abs_z < 0.90,
    prop_extreme < 0.10
  ) %>%
  arrange(mean_abs_z) %>%
  slice(1) %>%
  pull(subject)

if (length(normal_subject) == 0) {
  normal_subject <- subject_rank %>%
    filter(n_visits >= 3, prop_extreme < 0.15) %>%
    arrange(mean_abs_z) %>%
    slice(1) %>%
    pull(subject)
}

case_subjects <- c(normal_subject, abnormal_subject)

case_labels <- tibble(
  subject = case_subjects,
  case_group = c(
    "Representative low-deviation subject",
    "Highest-abnormality subject"
  )
)

dat_case <- dat_stan %>%
  filter(subject %in% case_subjects, region %in% key_regions) %>%
  left_join(case_labels, by = "subject") %>%
  mutate(
    region_clean = clean_region_name(region),
    case_group = factor(
      case_group,
      levels = c(
        "Highest-abnormality subject",
        "Representative low-deviation subject"
      )
    )
  )

p_case <- ggplot(dat_case, aes(x = age, y = z_C, group = region_clean)) +
  geom_line(aes(linetype = region_clean), linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ case_group, ncol = 1) +
  labs(
    subtitle = "Comparison of a representative low-deviation subject and the highest-abnormality subject",
    x = "Age at MRI visit",
    y = "Standardized deviation",
    linetype = "Region"
  ) +
  theme_pub(13) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(
  "manuscript/figures/Figure4_Subject_Case_Studies.pdf",
  p_case,
  width = 10.5,
  height = 8,
  bg = "white"
)

ggsave(
  "manuscript/figures/Figure4_Subject_Case_Studies.png",
  p_case,
  width = 10.5,
  height = 8,
  dpi = 500,
  bg = "white"
)

write_csv(
  subject_rank %>%
    filter(subject %in% case_subjects),
  "manuscript/tables/Table_Selected_Case_Study_Subjects.csv"
)

### SUPPLEMENTARY FIGURE S1: MEAN Z-SCORE CORTICAL ATLAS

p_brain_mean <- ggplot(atlas_df) +
  geom_brain(
    atlas = dk(),
    position = position_brain(hemi ~ view),
    aes(fill = mean_z)
  ) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-mean_lim, mean_lim),
    labels = scales::number_format(accuracy = 0.0001),
    na.value = "grey85"
  ) +
  labs(
    title = "Mean cortical standardized deviations",
    subtitle = "Posterior mean standardized deviations from the Bayesian spatial model",
    fill = "Mean z-score"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

p_region_mean_panel <- atlas_df %>%
  arrange(mean_z) %>%
  ggplot(aes(x = fct_reorder(region_hemi, mean_z), y = mean_z)) +
  geom_col(color = "black", width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = NULL,
    y = "Mean z-score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

p_s1_mean_atlas <- p_brain_mean + p_region_mean_panel +
  plot_layout(widths = c(2.2, 1))

ggsave(
  "manuscript/supplementary/figures/FigureS1_Mean_Cortical_Standardized_Deviations.pdf",
  p_s1_mean_atlas,
  width = 15,
  height = 7,
  bg = "white"
)

ggsave(
  "manuscript/supplementary/figures/FigureS1_Mean_Cortical_Standardized_Deviations.png",
  p_s1_mean_atlas,
  width = 15,
  height = 7,
  dpi = 600,
  bg = "white"
)

### SUPPLEMENTARY FIGURE S2: STANDARDIZED DEVIATION HISTOGRAMS

z_long <- dat_stan %>%
  select(subject, region, age, z_A, z_B, z_C) %>%
  pivot_longer(
    cols = c(z_A, z_B, z_C),
    names_to = "model",
    values_to = "z"
  ) %>%
  mutate(
    model = recode(
      model,
      z_A = "Independent cross-sectional",
      z_B = "Longitudinal non-spatial",
      z_C = "Bayesian subject-specific spatial"
    )
  )

p_s2_hist <- ggplot(z_long, aes(x = z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45, color = "black", alpha = 0.75) +
  stat_function(fun = dnorm, linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ model, ncol = 1) +
  labs(
    subtitle = "Dashed curve shows the standard normal density",
    x = "Standardized deviation",
    y = "Density"
  ) +
  theme_pub(13)

ggsave(
  "manuscript/supplementary/figures/FigureS2_Standardized_Deviation_Histograms.pdf",
  p_s2_hist,
  width = 8,
  height = 8,
  bg = "white"
)

ggsave(
  "manuscript/supplementary/figures/FigureS2_Standardized_Deviation_Histograms.png",
  p_s2_hist,
  width = 8,
  height = 8,
  dpi = 500,
  bg = "white"
)

### SUPPLEMENTARY FIGURE S3: QQ PLOT OF STANDARDIZED DEVIATIONS

p_s3_qq <- ggplot(z_long, aes(sample = z)) +
  stat_qq(alpha = 0.35, size = 0.8) +
  stat_qq_line(linewidth = 0.8) +
  facet_wrap(~ model, ncol = 1) +
  labs(
    x = "Theoretical quantiles",
    y = "Observed standardized deviations"
  ) +
  theme_pub(13)

ggsave(
  "manuscript/supplementary/figures/FigureS3_QQ_Plots_Standardized_Deviations.pdf",
  p_s3_qq,
  width = 8,
  height = 8,
  bg = "white"
)

ggsave(
  "manuscript/supplementary/figures/FigureS3_QQ_Plots_Standardized_Deviations.png",
  p_s3_qq,
  width = 8,
  height = 8,
  dpi = 500,
  bg = "white"
)

### SUPPLEMENTARY FIGURE S4: OBSERVED VS PREDICTED

obs_pred_long <- dat_stan %>%
  select(value_z, mu_A, mu_B, mu_C) %>%
  pivot_longer(
    cols = c(mu_A, mu_B, mu_C),
    names_to = "model",
    values_to = "predicted"
  ) %>%
  mutate(
    model = recode(
      model,
      mu_A = "Independent cross-sectional",
      mu_B = "Longitudinal non-spatial",
      mu_C = "Bayesian subject-specific spatial"
    )
  )

p_s4_obs_pred <- ggplot(obs_pred_long, aes(x = predicted, y = value_z)) +
  geom_point(alpha = 0.20, size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ model, ncol = 1) +
  labs(
    x = "Fitted value",
    y = "Observed standardized value"
  ) +
  theme_pub(13)

ggsave(
  "manuscript/supplementary/figures/FigureS4_Observed_vs_Predicted.pdf",
  p_s4_obs_pred,
  width = 8,
  height = 8,
  bg = "white"
)

ggsave(
  "manuscript/supplementary/figures/FigureS4_Observed_vs_Predicted.png",
  p_s4_obs_pred,
  width = 8,
  height = 8,
  dpi = 500,
  bg = "white"
)

### SUPPLEMENTARY FIGURE S5: MCMC TRACEPLOTS FOR KEY PARAMETERS

draws_key <- fit_C_stan$draws(
  variables = c("sigma", "sigma_b", "tau_u", "sigma_version")
)

p_s5_trace <- bayesplot::mcmc_trace(
  draws_key,
  pars = c("sigma", "sigma_b", "tau_u", "sigma_version")
) + theme_pub(13)

ggsave(
  "manuscript/supplementary/figures/FigureS5_MCMC_Traceplots_Key_Parameters.pdf",
  p_s5_trace,
  width = 10,
  height = 7,
  bg = "white"
)

ggsave(
  "manuscript/supplementary/figures/FigureS5_MCMC_Traceplots_Key_Parameters.png",
  p_s5_trace,
  width = 10,
  height = 7,
  dpi = 500,
  bg = "white"
)

### SUPPLEMENTARY FIGURE S6: SUBJECT ABNORMALITY BURDEN

p_s6_mean_abs <- ggplot(subject_dev_C, aes(x = mean_abs_z)) +
  geom_histogram(bins = 40, color = "black", alpha = 0.75) +
  labs(
    x = "Mean absolute standardized deviation",
    y = "Number of subjects"
  ) +
  theme_pub(13)

p_s6_prop_extreme <- ggplot(subject_dev_C, aes(x = prop_extreme)) +
  geom_histogram(bins = 40, color = "black", alpha = 0.75) +
  labs(
    x = "Proportion of observations with |z| > 1.96",
    y = "Number of subjects"
  ) +
  theme_pub(13)

p_s6_subject_burden <- p_s6_mean_abs / p_s6_prop_extreme

ggsave(
  "manuscript/supplementary/figures/FigureS6_Subject_Abnormality_Burden.pdf",
  p_s6_subject_burden,
  width = 8.5,
  height = 8,
  bg = "white"
)

ggsave(
  "manuscript/supplementary/figures/FigureS6_Subject_Abnormality_Burden.png",
  p_s6_subject_burden,
  width = 8.5,
  height = 8,
  dpi = 500,
  bg = "white"
)

### SUPPLEMENTARY CALIBRATION TABLE

write_csv(
  z_long %>%
    group_by(model) %>%
    summarise(
      mean_z = mean(z, na.rm = TRUE),
      var_z = var(z, na.rm = TRUE),
      tail_probability = mean(abs(z) > 1.96, na.rm = TRUE),
      .groups = "drop"
    ),
  "manuscript/supplementary/tables/TableS5_Calibration_Summary.csv"
)

