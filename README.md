# Bayesian Longitudinal Spatial Normative Modeling

This repository contains the code, simulation framework, and manuscript materials for a Bayesian longitudinal spatial normative modeling framework for individualized structural neuroimaging analysis.

The proposed methodology integrates:

- longitudinal repeated-measure dependence,
- subject-specific latent spatial deviation processes,
- and region-level Bayesian normative modeling

within a unified hierarchical framework for individualized deviation estimation from structural MRI data.

The repository includes:

- simulation studies across multiple longitudinal and spatial scenarios,
- OASIS-3 structural MRI application workflows,
- Stan model implementations,
- manuscript and supplementary materials,
- and reproducible figure/table generation pipelines.

---

# Repository Structure

```text
R/                         reusable utility functions
stan/                      Stan model files

simulations/
  scripts/                 simulation workflows
  figures/                 simulation figures
  results/                 simulation summary tables

oasis_application/
  preprocessing/           OASIS-3 preprocessing workflows
  modeling/                model fitting and posterior extraction
  outputs/                 figure and table generation scripts

manuscript/
  figures/                 main manuscript figures
  supplementary/           supplementary figures and tables


Main Methodological Components

The framework jointly models:

region-specific normative trajectories,
subject-level repeated measurements,
and spatially structured individualized deviation maps.

The Bayesian model includes:

region-level fixed effects,
subject-specific random intercepts,
and latent spatial deviation processes with structured covariance.

Simulation studies evaluate:

deviation-map recovery,
calibration,
estimation accuracy,
and abnormality detection

under multiple longitudinal and spatial data-generating scenarios.


Data Availability

The OASIS-3 neuroimaging data used in this work are publicly available through the OASIS project under their data use agreement:

https://www.oasis-brains.org/

Raw participant-level data are not distributed in this repository.



Software Requirements

Main R packages used include:

cmdstanr
posterior
dplyr
tidyr
ggplot2
lme4
patchwork
ggseg
bayesplot
viridis

Stan models were fit using CmdStan.

Reproducibility

All scripts are intended to be run from the repository root directory.

Typical workflow:

Rscript oasis_application/preprocessing/01_prepare_oasis_data.R

Rscript oasis_application/modeling/03_build_modeling_data_and_adjacency.R

Rscript oasis_application/modeling/04_fit_oasis_models.R

Rscript oasis_application/modeling/05_extract_oasis_deviation_scores.R

Rscript oasis_application/outputs/06_make_oasis_tables.R

Rscript oasis_application/outputs/07_make_oasis_figures.R


Simulation summaries can be regenerated using:

Rscript simulations/scripts/summarize_simulation_results.R


Manuscript Status

Methodological manuscript and supplementary materials currently included in repository.

Additional package development and methodological extensions are ongoing.
