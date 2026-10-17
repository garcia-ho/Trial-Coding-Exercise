# Phase I dose-finding practice

A small R project for practicing phase I dose finding with the continual
reassessment method (CRM), BMA-CRM, and two-drug combination BOIN.

## Files

- `CRM/CRM.Rmd`: CRM, skeleton, MCMC, and BMA-CRM experiment
- `CRM/R/crm_functions.R`: commented CRM functions
- `Drug_combo_Boin/Drug_combo_Boin.Rmd`: combination-BOIN simulation study
- `Drug_combo_Boin/R/drug_combo_boin_functions.R`: commented combination-BOIN
  functions
- `Drug_combo_Boin/Drug_combo_Boin.html`: rendered combination-BOIN report

Open `phase1-dose-finding.Rproj`, then open either R Markdown file. Run an
individual chunk with its green play button or render the full notebook with
**Knit**.

From the project root, the combination report can also be rebuilt with:

```r
rmarkdown::render("Drug_combo_Boin/Drug_combo_Boin.Rmd")
```

The exercises use base R plus `knitr` and `rmarkdown`; no dose-finding package
is required. The implementations are educational and are not validated
software for conducting a clinical trial.
