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
- `Lyu_2018/simulation_exploration.Rmd`: saved-result allocation and scenario analysis, with a conceptual copula stress-test proposal
- `Lyu_2018/simulation_exploration.html`: rendered exploration report
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
rmarkdown::render("Lyu_2018/simulation_exploration.Rmd")
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
Go decisions, all 15 exits, and the week-56 finish; all 42 checks pass.
The report documents three source limitations:

- At week 40, cycle 3 D2, Table 2 prints 1/6; its patient rows give 1/3.
- Figure 2's 2-DLT/2-patient decision differs from the appendix's unqualified
  posterior safety rule. `safety_min_n = 3` reproduces the printed figure;
  either convention gives the same Table 2 dose assignments.
- PAVA weights and nearest-target ties within an ordering are incompletely
  specified. The MTS function requires these choices explicitly and the
  report shows their effect. Table 2 does not publish a final MTS.

The Table 2 section is a deterministic replay of the supplied patient
outcomes. A Table 3 simulation comparison is appended in the same Rmd.
See Lyu, Curran and Ji (2018),
[DOI 10.1200/PO.18.00020](https://doi.org/10.1200/PO.18.00020), Table 2
(pp. 6-8) and Appendix (pp. 16-18).

## Dependencies

The exercises use base R plus `knitr` and `rmarkdown`; no dose-finding package
is required. Rendering also needs Pandoc, included with RStudio or Quarto.
The implementations are educational and are not validated
software for conducting a clinical trial.

## Table 3: one simulation file, one data file, the same Rmd

Use a **local VS Code window and its local terminal**, or macOS Terminal.
Keep the Mac as the editing copy; Unity runs the larger simulation. The existing
Table 2 Rmd now includes Table 3 after its original sections.

- `Lyu_2018/simulate.R`: all Table 3 simulation and summary functions. It calls
  the existing `run_basyc()` and `select_mts()` functions.
- `Lyu_2018/data/table3.rds`: one saved dataset, including settings and raw records.
- `Lyu_2018/reproduce_table2.Rmd`: Table 2 followed by the saved Table 3 comparisons.
- `unity.sh`: upload, submit the parallel job, download the RDS, then render locally.

From a **local** terminal:

```sh
cd "/Users/jesush/Documents/OSU/Clinical Trial/Trial Coding Exercise"
bash unity.sh 10000 auto 2018
```

Arguments are trials **per scenario**, CPU cores (`auto` or a number), and seed.
The run covers all six scenarios: 10,000 each gives 60,000 trials total.
`auto` checks healthy nodes in the `stat` partition available to your account and requests up to
32 currently idle cores on one node; Slurm makes the actual allocation.
For a quick server check use `bash unity.sh 100 auto 2018`. Enter password/Duo privately. Everything after
login runs automatically. `simulate_table3()` uses `parallel::mclapply()` over trial batches, so six
scenarios can use more than six cores. R never exceeds its Slurm allocation.
The RDS and report include total wall-clock simulation runtime, excluding the
scheduler wait and upload/download time.

The remote working folder is `~/trial-coding-exercise`. Files persist there.
Rsync updates changed source files and skips generated results, the saved RDS,
and installed Mac packages. There are no per-run folders.
Each successful new simulation replaces `data/table3.rds` and updates the same
HTML report. The job's single diagnostic log stays on Unity as
`Lyu_2018/table3.log`; it is not downloaded.

Keep the terminal open until download/rendering finishes. Do not edit and
upload the folder again while its job is running. If the connection drops,
the submitted job can continue; check `squeue -u he.2774` after logging in.
After confirming successful completion, retrieve the existing data without
rerunning the simulation:

```sh
scp asc-unity:trial-coding-exercise/Lyu_2018/data/table3.rds Lyu_2018/data/table3.rds
Rscript --vanilla -e 'rmarkdown::render("Lyu_2018/reproduce_table2.Rmd")'
```

Slurm defaults are partition `stat`, eight GB RAM and one hour; the R module
is `r/4.5.1-5zezbzn`. Change the clearly marked values in `unity.sh` if needed.
For example, prefix a command with `UNITY_TIME=00:30:00`. Simulations run on
compute nodes, as required by the [Unity walk-through](https://osu.teamdynamix.com/TDClient/1929/ASC/KB/Article/61538/Unity-Walk-Through).

The simulation default requires strict **n > 3** for intratrial safety,
implemented by `safety_min_n = 4`. The same minimum governs DCR, global
termination and permanent dose exclusion/DU. It is an explicit simulation
hypothesis; Figure 2's printed 3/3 decision uses a different minimum.
Final MTS posterior screening has a separate `mts_screen_min_n = 0`
argument, retaining the appendix's unqualified screen. Both settings are saved
with new results. PAVA weights remain equal and order enumeration is
unchanged. `nearest_tie = "higher"` uses the original higher-position rule
for nearest-target ties, including in the existing fallback. The
target-dependent convention is not used.

`global_monitor = "after_first_sequence"` separates global safety from DCR:
DCR operates from week 8, while global termination and permanent exclusion
both begin when cohort 1 completes all three cycles (week 24). Outcomes due
are observed first; active global safety checks precede further treatment
and enrollment. A final visit after all treatment is complete does not count
as early termination, but its data still inform final safety and MTS selection.
The RDS records global activation and the first stopping cell's counts/risk.
This is an explicit sensitivity experiment; the paper describes safety
monitoring throughout the trial. Delaying exclusion with termination avoids
an incoherent state with no admissible dose but an instruction to continue.

`no_mts_fallback = "restricted_revote"` handles the case in which every
original winner is increasing: retain the fitted TPCs and repeat the
nearest-target vote among observed, safe, nonincreasing sequences. The final
conservative frequency tie rule remains unchanged. Early stopping and an
empty safe nonincreasing set still give no MTS. This fallback is an explicit
completion of an appendix edge case, not a verified author implementation.

Knitting the Rmd **only reads saved Table 3 data**. It lists the functions
called and compares selection, completed-sequence DLT/allocation, dose/cycle
counts and early termination with transcribed Table 3 values. No comparison
with earlier simulation runs is displayed. Sections 5.2 and 6.3 state the
NO-MTS rule and the saved data record each trial's fallback use. The report
also gives runtime, settings and remaining discrepancies with the paper.
The Table 2 allocation replay is unchanged; the current illustration
selects 221 without using a fallback. Table 2 publishes no MTS.
The saved RDS and report give the current run's settings, runtime, server
job details and comparison with the paper for all six scenarios. These
details are read from the saved run rather than duplicated here.

For direct access next time: `ssh asc-unity`, or in VS Code use
**Cmd+Shift+P -> Remote-SSH: Connect to Host... -> asc-unity**. A remote VS Code
terminal runs on Unity; the `unity.sh` command above belongs in a **local**
terminal. For the present upload-and-render workflow, a local window is simpler.

To let Codex reuse a connection that you authenticate privately in a local
Terminal, start a shared SSH connection and leave that Terminal open:

```sh
ssh -M -S ~/.ssh/asc-unity.sock -o ServerAliveInterval=30 -N asc-unity
```

Enter password/Duo in Terminal, then tell Codex that it is connected. A blank
Terminal after authentication is normal. No password is shared or saved.
The workflow can reuse that connection with
`UNITY_SOCKET=/Users/jesush/.ssh/asc-unity.sock bash unity.sh 10000 auto 2018`.
Press Ctrl+C in the authentication Terminal when the work is finished.

## Server terminal

zsh 5.5.1 is already installed on Unity at `/usr/bin/zsh`. New interactive SSH
sessions and remote VS Code terminals are configured to use it. The original
`.bashrc` was backed up on Unity; `.zshrc` was preserved. Slurm and noninteractive
SSH commands continue using bash. The managed account login shell is unchanged;
ASCTech must change that directory-service setting if desired.

## Saved-result exploration

`Lyu_2018/simulation_exploration.Rmd` reads the existing `data/table3.rds`
without sourcing the simulator or generating trial outcomes. It displays
recommendation categories, dose/cycle allocation bars and heatmaps, completed
patient transitions, and true-MTS coverage versus selection. Scenario-specific
interpretations distinguish exploration, selection and patient exposure.
The final section proposes one- and two-parameter copula stress tests; no
copula generator or new simulation is implemented. Knit this report separately
when the saved RDS changes; `unity.sh` continues to render the reproduction report.
