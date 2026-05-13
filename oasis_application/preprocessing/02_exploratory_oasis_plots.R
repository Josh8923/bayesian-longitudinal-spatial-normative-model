#---------------------------------------- ========================
######### OASIS-3 Exploratory Plots
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)

source("R/utilities.R")

ensure_dir("oasis_application/outputs/qc_figures")

dat_long <- readRDS("oasis_application/outputs/oasis3_model_ready_long.rds")
dat_wide <- readRDS("oasis_application/outputs/oasis3_model_ready_wide.rds")
dat_long_2plus <- readRDS("oasis_application/outputs/oasis3_model_ready_long_2plus.rds")
dat_wide_2plus <- readRDS("oasis_application/outputs/oasis3_model_ready_wide_2plus.rds")

subject_summary <- read_csv(
  "oasis_application/outputs/oasis3_subject_followup_summary.csv",
  show_col_types = FALSE
)

### MRI VISITS PER SUBJECT

p_visits <- subject_summary %>%
  ggplot(aes(x = n_visits)) +
  geom_bar(color = "black") +
  labs(
    title = "Distribution of MRI visits per subject",
    x = "Number of MRI visits",
    y = "Number of subjects"
  ) +
  theme_pub(12)

print(p_visits)

ggsave(
  "oasis_application/outputs/qc_figures/plot_visits_per_subject.png",
  p_visits,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)

### AGE DISTRIBUTION

p_age <- dat_wide %>%
  ggplot(aes(x = age)) +
  geom_histogram(bins = 30, color = "black") +
  labs(
    title = "Age distribution across MRI sessions",
    x = "Age at visit",
    y = "Number of MRI sessions"
  ) +
  theme_pub(12)

print(p_age)

ggsave(
  "oasis_application/outputs/qc_figures/plot_age_distribution.png",
  p_age,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)

### OBSERVED LONGITUDINAL TRAJECTORIES

selected_plot_regions <- c(
  "left_hippocampus_volume",
  "right_hippocampus_volume",
  "left_amygdala_volume",
  "right_amygdala_volume",
  "lh_entorhinal_thickness",
  "rh_entorhinal_thickness",
  "lh_parahippocampal_thickness",
  "rh_parahippocampal_thickness"
)

selected_plot_regions <- selected_plot_regions[
  selected_plot_regions %in% levels(dat_long$region)
]

p_traj <- dat_long_2plus %>%
  filter(region %in% selected_plot_regions) %>%
  ggplot(aes(x = age, y = value, group = subject)) +
  geom_line(alpha = 0.12) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ region, scales = "free_y", ncol = 2) +
  labs(
    title = "Observed longitudinal trajectories for selected structural regions",
    x = "Age at visit",
    y = "ICV-normalized volume or cortical thickness"
  ) +
  theme_pub(12)

print(p_traj)

ggsave(
  "oasis_application/outputs/qc_figures/plot_selected_region_trajectories.png",
  p_traj,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

### EMPIRICAL REGION CORRELATION HEATMAP

region_wide_only <- dat_wide_2plus %>%
  select(all_of(unique(as.character(dat_long_2plus$region)))) %>%
  select(where(is.numeric))

region_cor <- cor(region_wide_only, use = "pairwise.complete.obs")

write.csv(
  region_cor,
  "oasis_application/outputs/qc_figures/oasis3_region_correlation_matrix.csv"
)

region_cor_long <- as.data.frame(region_cor) %>%
  mutate(region1 = rownames(.)) %>%
  pivot_longer(
    cols = -region1,
    names_to = "region2",
    values_to = "correlation"
  )

p_cor <- ggplot(region_cor_long, aes(x = region1, y = region2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1)) +
  labs(
    title = "Empirical correlation among selected brain regions",
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

print(p_cor)

ggsave(
  "oasis_application/outputs/qc_figures/plot_region_correlation_heatmap.png",
  p_cor,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

