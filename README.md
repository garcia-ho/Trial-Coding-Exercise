# Phase I dose-finding practice

A small R project for practicing the continual reassessment method (CRM),
indifference-interval skeleton calibration, MCMC, and BMA-CRM.

## Files

- `CRM.Rmd`: executable CRM tutorial and simulation study
- `R/crm_functions.R`: reusable, commented functions called by the notebook
- `data/`: exercise datasets

Open `phase1-dose-finding.Rproj`, then open `CRM.Rmd`. Run an individual
chunk with its green play button or render the full notebook with **Knit**.

The exercise uses base R plus `knitr` and `rmarkdown`; no dose-finding package
is required. The implementation is educational and is not validated software
for conducting a clinical trial.
