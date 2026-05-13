
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(readr)

# Load simulation results
res1 <- readRDS("simulations/results/sim_results_Scenario_1_NoSpatial.rds")
res2 <- readRDS("simulations/results/sim_results_Scenario_2_ModerateSpatial.rds")
res3 <- readRDS("simulations/results/sim_results_Scenario_3_StrongSpatial.rds")
res4 <- readRDS("simulations/results/sim_results_Scenario_4_VariableVisits.rds")
res5 <- readRDS("simulations/results/sim_results_Scenario_5_MissingFollowup.rds")
res6 <- readRDS("simulations/results/sim_results_Scenario_6_NonlinearMean.rds")

all_res <- list(res1, res2, res3, res4, res5, res6)

scenario_labels <- c(
  "Scenario_1_NoSpatial" = "No spatial dependence",
  "Scenario_2_ModerateSpatial" = "Moderate spatial dependence",
  "Scenario_3_StrongSpatial" = "Strong spatial dependence",
  "Scenario_4_VariableVisits" = "Variable visit schedule",
  "Scenario_5_MissingFollowup" = "Missing follow-up",
  "Scenario_6_NonlinearMean" = "Nonlinear mean structure",
  "Scenario_6_Nonlinear" = "Nonlinear mean structure"
)

scenario_order <- c(
  "No spatial dependence",
  "Moderate spatial dependence",
  "Strong spatial dependence",
  "Variable visit schedule",
  "Missing follow-up",
  "Nonlinear mean structure"
)

model_labels <- c(
  "Model_A_Independent_CS" = "Independent cross-sectional model",
  "Model_B_Longitudinal_NoSpatial" = "Longitudinal non-spatial model",
  "Model_C_Bayesian_Longitudinal_Spatial" = "Bayesian longitudinal spatial model"
)

model_order <- c(
  "Independent cross-sectional model",
  "Longitudinal non-spatial model",
  "Bayesian longitudinal spatial model"
)

# Master summary
master_summary <- bind_rows(lapply(all_res, function(x) x$summary)) %>%
  mutate(
    `Simulation scenario` = recode(scenario, !!!scenario_labels),
    `Model` = recode(model, !!!model_labels),
    `Simulation scenario` = factor(`Simulation scenario`, levels = scenario_order),
    `Model` = factor(`Model`, levels = model_order)
  ) %>%
  arrange(`Simulation scenario`, `Model`) %>%
  transmute(
    `Simulation scenario`,
    `Model`,
    `Mean bias` = mean_bias_mean,
    `Mean squared error` = mean_mse_mean,
    `Deviation map MSE` = map_mse_mean,
    `z-score mean` = z_mean_mean,
    `z-score variance` = z_var_mean,
    `Tail probability` = z_tail_prop_mean
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# ## Table 1: Estimation accuracy
table_accuracy <- master_summary %>%
  select(
    `Simulation scenario`,
    `Model`,
    `Mean bias`,
    `Mean squared error`,
    `Deviation map MSE`
  )

table_accuracy %>%
  kbl(
    caption = "Table 1. Estimation accuracy across simulation scenarios.
    Lower mean squared error and deviation map MSE indicate better performance.",
    booktabs = TRUE,
    longtable = TRUE,
    align = "llrrr"
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    full_width = FALSE,
    font_size = 9
  ) %>%
  column_spec(1, bold = TRUE, width = "3.6cm") %>%
  column_spec(2, width = "4.5cm") %>%
  collapse_rows(columns = 1, valign = "top")

# ##Table 2: Calibration
table_calibration <- master_summary %>%
  select(
    `Simulation scenario`,
    `Model`,
    `z-score mean`,
    `z-score variance`,
    `Tail probability`
  )

table_calibration %>%
  kbl(
    caption = "Table 2. Calibration performance across simulation scenarios.
    Target values are 0 for z-score mean, 1 for z-score variance,
    and 0.05 for tail probability.",
    booktabs = TRUE,
    longtable = TRUE,
    align = "llrrr"
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    full_width = FALSE,
    font_size = 9
  ) %>%
  column_spec(1, bold = TRUE, width = "3.6cm") %>%
  column_spec(2, width = "4.5cm") %>%
  collapse_rows(columns = 1, valign = "top")

###  Raw data for figures
raw_all <- bind_rows(lapply(all_res, function(x) x$raw)) %>%
  mutate(
    `Simulation scenario` = recode(scenario, !!!scenario_labels),
    `Model` = recode(model, !!!model_labels),
    `Simulation scenario` = factor(`Simulation scenario`, levels = scenario_order),
    `Model` = factor(`Model`, levels = model_order)
  )

model_colors <- c(
  "Independent cross-sectional model" = "#4D4D4D",
  "Longitudinal non-spatial model" = "#0072B2",
  "Bayesian longitudinal spatial model" = "#D55E00"
)

publication_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  )

# Figure 1: Estimation accuracy
p_map <- ggplot(raw_all, aes(x = `Model`, y = map_mse, fill = `Model`)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, linewidth = 0.3) +
  facet_wrap(~ `Simulation scenario`, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Deviation map reconstruction accuracy",
    x = NULL,
    y = "Deviation map mean squared error"
  ) +
  publication_theme

p_mse <- ggplot(raw_all, aes(x = `Model`, y = mean_mse, fill = `Model`)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, linewidth = 0.3) +
  facet_wrap(~ `Simulation scenario`, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Global estimation accuracy",
    x = NULL,
    y = "Mean squared error"
  ) +
  publication_theme

fig_accuracy <- p_map / p_mse +
  plot_annotation(
    #title = "Figure 1. Estimation accuracy and deviation map reconstruction across simulation scenarios",
    caption = "Boxplots summarize Monte Carlo distributions across simulation replicates. Lower values indicate better performance."
  )

print(fig_accuracy)

ggsave(
  "figure1_accuracy_across_scenarios.png",
  fig_accuracy,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

# Figure 2: Calibration
p_zmean <- ggplot(raw_all, aes(x = `Model`, y = z_mean, fill = `Model`)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ `Simulation scenario`, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "z-score centering",
    x = NULL,
    y = "z-score mean"
  ) +
  publication_theme

p_zvar <- ggplot(raw_all, aes(x = `Model`, y = z_var, fill = `Model`)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ `Simulation scenario`, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "z-score dispersion",
    x = NULL,
    y = "z-score variance"
  ) +
  publication_theme

p_tail <- ggplot(raw_all, aes(x = `Model`, y = z_tail_prop, fill = `Model`)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, linewidth = 0.3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ `Simulation scenario`, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Tail calibration",
    x = NULL,
    y = "P(|z| > 1.96)"
  ) +
  publication_theme

fig_calibration <- (p_zmean / p_zvar / p_tail) +
  plot_annotation(
   # title = "Figure 2. Calibration under varying spatial and longitudinal data-generating mechanisms",
    caption = "Dashed reference lines indicate target values: 0 for z-score mean, 1 for z-score variance, and 0.05 for tail probability."
  )

print(fig_calibration)

ggsave(
  "figure2_calibration_across_scenarios.png",
  fig_calibration,
  width = 14,
  height = 14,
  dpi = 300,
  bg = "white"
)

# save clean CSV tables
write_csv(table_accuracy, "table1_estimation_accuracy.csv")
write_csv(table_calibration, "table2_calibration.csv")



ggsave(
  "simulations/figures/figure2_calibration_across_scenarios.pdf",
  fig_calibration,
  width = 14,
  height = 14,
  bg = "white"
)

ggsave(
  "simulations/figures/figure1_accuracy_across_scenarios.pdf",
  fig_accuracy,
  width = 14,
  height = 10,
  bg = "white"
)
