# Phase I dose-finding practice

A shared R project for practicing phase I dose finding with the continual
reassessment method (CRM), BMA-CRM, two-drug combination BOIN, and the
multicycle BaSyc design of Lyu et al. (2018).

## Files

- `CRM/CRM.Rmd`: CRM, skeleton, MCMC, and BMA-CRM experiment
- `CRM/R/crm_functions.R`: commented CRM functions
- `Drug_combo_Boin/Drug_combo_Boin.Rmd`: combination-BOIN simulation study
- `Drug_combo_Boin/R/drug_combo_boin_functions.R`: commented combination-BOIN
  functions
- `Drug_combo_Boin/Drug_combo_Boin.html`: rendered combination-BOIN report
- `Lyu_2018/reproduce_table2.Rmd`: BaSyc formulas, Table 2 replay and checks
- `Lyu_2018/reproduce_table2.html`: rendered BaSyc report
- `Lyu_2018/R/`: commented mTPI-2, allocation, safety, MTS and table functions
- `Lyu_2018/validate.R`: independently runnable reproduction checks
- `Lyu_2018/data/`: transcribed outcomes, expected doses and printed count tables
- `Lyu_2018/results/`: generated assignments, count comparisons and validation tables

## One project, independent exercises

Open `phase1-dose-finding.Rproj`, then open an exercise's R Markdown file. Run an
individual chunk with its green play button or render the full notebook with
**Knit**. The root project and its `renv` environment are shared; exercise
folders do not need their own `.Rproj` files or render scripts.

From the project root, rebuild individual reports directly:

```r
rmarkdown::render("CRM/CRM.Rmd")
rmarkdown::render("Drug_combo_Boin/Drug_combo_Boin.Rmd")
rmarkdown::render("Lyu_2018/reproduce_table2.Rmd")
```

## Lyu 2018: run and review

The Lyu exercise also runs directly from its own folder. Its paths support
both the shared project root and the exercise folder, without changing the
working directory or relying on another exercise's functions.

From the shared project root:

```r
source("Lyu_2018/validate.R")
```

Or from a terminal at the shared project root:

```sh
Rscript --vanilla Lyu_2018/validate.R
```

If the working directory is already `Lyu_2018`:

```r
source("validate.R")
rmarkdown::render("reproduce_table2.Rmd")
```

The functions in `R/` contain the statistical logic. The R Markdown file
shows how to call them and displays the results. The validation script
recreates the CSVs in `Lyu_2018/results/`; it is retained so checks can run
without rendering the report. Inputs in `data/` are kept separate from
generated outputs in `results/`.

The reproduction matches all 45 Table 2 dose assignments, all 30 continuation
Go decisions, all 15 exits, and the week-56 finish; all 21 checks pass.
The report documents three source limitations:

- At week 40, cycle 3 D2, Table 2 prints 1/6; its patient rows give 1/3.
- Figure 2's 2-DLT/2-patient decision differs from the appendix's unqualified
  posterior safety rule. `safety_min_n = 3` reproduces the printed figure;
  either convention gives the same Table 2 dose assignments.
- PAVA weights and nearest-target ties within an ordering are incompletely
  specified. The MTS function requires these choices explicitly and the
  report shows their effect. Table 2 does not publish a final MTS.

This is a deterministic replay of the supplied patient outcomes, not a
reproduction of the paper's operating-characteristic simulations.
See Lyu, Curran and Ji (2018),
[DOI 10.1200/PO.18.00020](https://doi.org/10.1200/PO.18.00020), Table 2
(pp. 6-8) and Appendix (pp. 16-18).

## Dependencies

The exercises use base R plus `knitr` and `rmarkdown`; no dose-finding package
is required. Rendering also needs Pandoc, included with RStudio or Quarto.
The implementations are educational and are not validated
software for conducting a clinical trial.
