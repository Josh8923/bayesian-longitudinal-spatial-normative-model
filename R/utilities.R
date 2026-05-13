#---------------------------------------- ========================
######### Shared utilities
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(stringr)
library(ggplot2)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 2),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = base_size - 3),
      legend.position = "bottom"
    )
}

clean_region_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_replace_all("^lh ", "Left ") %>%
    str_replace_all("^rh ", "Right ") %>%
    str_replace_all("^left ", "Left ") %>%
    str_replace_all("^right ", "Right ") %>%
    str_replace_all("inferiortemporal", "inferior temporal") %>%
    str_replace_all("middletemporal", "middle temporal") %>%
    str_replace_all("posteriorcingulate", "posterior cingulate") %>%
    str_replace_all("temporalpole", "temporal pole") %>%
    str_replace_all("parahippocampal", "parahippocampal") %>%
    str_replace_all("entorhinal", "entorhinal") %>%
    str_to_sentence()
}

get_base_region <- function(x) {
  x %>%
    str_replace("^lh_", "") %>%
    str_replace("^rh_", "") %>%
    str_replace("^left_", "") %>%
    str_replace("^right_", "") %>%
    str_replace("_volume$", "") %>%
    str_replace("_thickness$", "")
}

get_mean_vector <- function(draws, prefix, n) {
  cols <- paste0(prefix, "[", seq_len(n), "]")
  colMeans(draws[, cols, drop = FALSE])
}
