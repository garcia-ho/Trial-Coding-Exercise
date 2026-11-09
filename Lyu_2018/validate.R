# From the shared project root: Rscript --vanilla Lyu_2018/validate.R
# From the exercise folder:     Rscript --vanilla validate.R
# Tests use base R; the expected doses/totals are NOT inputs to run_basyc().
lyu_dir <- if (file.exists("Lyu_2018/R/mtpi2.R")) "Lyu_2018" else "."
for (file in c("mtpi2", "basyc", "mts", "table2_helpers")) {
  source(file.path(lyu_dir, "R", paste0(file, ".R")))
}

outcomes <- read_table2_outcomes()
fit <- run_basyc(outcomes)
doses <- compare_table2_doses(fit)
totals <- compare_table2_totals(fit)
checks <- data.frame(check = character(), passed = logical())
check <- function(name, condition) {
  ok <- isTRUE(condition)
  checks[nrow(checks) + 1L, ] <<- list(name, ok)
  if (!ok) stop("FAILED: ", name)
}

check("All 45 independently transcribed Table 2 dose assignments match",
      nrow(doses) == 45 && all(doses$match))
check("All 30 continuation decisions are Go; all 15 final decisions are Exit",
      sum(fit$continuation$decision == "G") == 30 &&
      sum(fit$continuation$decision == "Exit") == 15)
check("Every cycle ends at the Table 2 calendar week",
      all(fit$observed$end_week ==
            8 * (fit$observed$cohort + fit$observed$cycle - 1)))
check("The replay ends at week 56 with all 45 outcomes observed",
      fit$end_week == 56 && nrow(fit$observed) == 45 && !fit$stopped)
check("All nine final toxicity numerators and denominators match",
      all(subset(totals, week == 56)$match))
diff <- subset(totals, !match)
check("The sole printed-total discrepancy is week 40, C3 D2: 1/6 versus 1/3",
      nrow(diff) == 1 && diff$week == 40 && diff$cycle == 3 && diff$dose == 2 &&
      diff$printed == "1/6" && diff$computed == "1/3")

check("mTPI-2 gives E at 0/3, S at 1/3, D at 3/6",
      identical(vapply(c(0, 1, 3), function(x) {
        mtpi2_decision(x, if (x == 3) 6 else 3)
      }, character(1)), c("E", "S", "D")))
ref <- figure2_reference()
appendix_table <- decision_tables()
figure_table <- decision_tables(settings = basyc_settings(safety_min_n = 3))
bad <- which(ref$DAR != appendix_table$DAR | ref$DCR != appendix_table$DCR)
check("Appendix versus Figure 2 differs only at 2 DLTs / 2 patients",
      length(bad) == 1 && ref$n[bad] == 2 && ref$x[bad] == 2)
check("The explicit minimum-n=3 convention matches all 135 cells of both Figure 2 panels",
      all(ref$DAR == figure_table$DAR) && all(ref$DCR == figure_table$DCR))
figure_fit <- run_basyc(outcomes, settings = basyc_settings(safety_min_n = 3))
check("The Figure 2 convention also reproduces all 45 Table 2 doses",
      all(compare_table2_doses(figure_fit)$match))
check("n>=3 preserves every Table 2 assignment, observed outcome and continuation decision",
      identical(figure_fit$assignments, fit$assignments) &&
      identical(figure_fit$observed, fit$observed) &&
      identical(figure_fit$continuation, fit$continuation) &&
      figure_fit$end_week == fit$end_week)

# A future outcome must not alter an assignment already made at week 16.
future_changed <- outcomes
for (id in seq_len(nrow(outcomes))) for (j in seq_len(ncol(outcomes))) {
  if (8 * (ceiling(id / 3) + j - 1) > 16) future_changed[id, j] <- 1L
}
changed_fit <- run_basyc(future_changed)
check("No future outcomes leak into decisions made by week 16",
      identical(subset(fit$assignments, week <= 16),
                subset(changed_fit$assignments, week <= 16)))
all_toxic <- run_basyc(matrix(1L, 3, 3))
check("Three DLTs at the lowest dose stop the trial at week 8",
      all_toxic$stopped && all_toxic$end_week == 8 &&
      nrow(all_toxic$assignments) == 3)
zero_toxic <- run_basyc(matrix(0L, 15, 3), settings = basyc_settings(n_dose = 3))
check("No-DLT replay respects dose boundaries and the first-cohort default",
      all(zero_toxic$assignments$dose %in% 1:3) &&
      all(subset(zero_toxic$assignments, cohort == 1)$dose == 1))

# Test permanent exclusion even if later counts would reduce the tail risk.
counts <- list(x = matrix(0L, 3, 3), n = matrix(0L, 3, 3))
counts$x[2, 1] <- counts$n[2, 1] <- 3L
s <- update_safety(counts, matrix(FALSE, 3, 3), basyc_settings())
check("Unsafe D2 in C1 excludes D2-D3 in all later cycles",
      all(s$excluded[2:3, ]) && !any(s$excluded[1, ]) && !s$stop_trial)
counts$n[2, 1] <- 30L
s2 <- update_safety(counts, s$excluded, basyc_settings())
check("Exclusion remains permanent after later posterior risk decreases",
      identical(s$excluded, s2$excluded))
check("DCR checks the lowest NEXT-cycle cell",
      dcr(1, list(x = matrix(c(0, 0, 0, 3, 0, 0), 3, 2),
                  n = matrix(c(0, 0, 0, 3, 0, 0), 3, 2)), basyc_settings()) == "NG")

check("Equal-weight PAVA pools an order violation",
      isTRUE(all.equal(pava(c(0.1, 0.3, 0.2, 0.4)), c(0.1, 0.25, 0.25, 0.4))))
check("Weighted PAVA uses the block's total weight",
      isTRUE(all.equal(pava(c(0.1, 0.3, 0.2, 0.4), c(1, 1, 3, 1)),
                       c(0.1, 0.225, 0.225, 0.4))))
check("PAVA rechecks larger blocks after a merge",
      isTRUE(all.equal(pava(c(0.5, 0.4, 0.1)), rep(1/3, 3))))
lower <- select_mts(fit, weights = "equal", nearest_tie = "lower")
higher <- select_mts(fit, weights = "equal", nearest_tie = "higher")
check("MTS considers 2! x 3! = 12 possible intensity-respecting orders",
      nrow(lower$per_order) == 12 && sum(lower$votes$frequency) == 12)
check("Explicit within-order tie conventions expose different MTS outputs",
      lower$recommendation == "2,1,1" && higher$recommendation == "2,2,1")
revised_mts <- select_mts(figure_fit, weights = "equal", nearest_tie = "higher",
                          screen_min_n = 0)
check("The restored PAVA rule retains Table 2 MTS 221 without using the fallback",
      revised_mts$recommendation == higher$recommendation &&
      revised_mts$votes$frequency[revised_mts$votes$sequence == "2,2,1"] == 6 &&
      !revised_mts$fallback_used)
for (r in seq_len(nrow(lower$per_order))) {
  ord <- strsplit(lower$per_order$ordering[r], " < ", fixed = TRUE)[[1]]
  stopifnot(all(diff(lower$fitted[r, ord]) >= -1e-12))
}

# A small final-selection fixture, not a simulated or allocated trial:
# 111 has 0/9, 211 has 3/9, 121 has 1/9. All dose/cycle cells pass screening.
# In order 211<121, PAVA pools their TPCs; the restored rule selects the
# higher position within that ordering, rather than imposing a clinical tie.
fixture_paths <- rbind(matrix(1L, 3, 3), matrix(rep(c(2L,1L,1L), each=3), 3),
                       matrix(rep(c(1L,2L,1L), each=3), 3))
fixture_y <- matrix(0L, 9, 3)
fixture_y[cbind(4:6, 1:3)] <- 1L
fixture_y[7,2] <- 1L
fixture_observed <- expand.grid(patient=1:9, cycle=1:3)
fixture_observed$dose <- as.vector(fixture_paths)
fixture_observed$tox <- as.vector(fixture_y)
fixture <- list(observed=fixture_observed, n_patient=9L, n_cycle=3L,
                settings=basyc_settings(safety_min_n=3), stopped=FALSE,
                counts=dose_cycle_counts(fixture_observed,3,3),
                snapshots=list(list(excluded=matrix(FALSE,3,3))))
legacy <- select_mts(fixture, weights="equal", nearest_tie="higher")
check("Original within-order PAVA ties give one vote each to 211 and 121",
      nrow(legacy$per_order)==2 &&
      identical(sort(legacy$per_order$selected), c("1,2,1","2,1,1")))

# A deterministic missing-winner example, not a Monte Carlo trial:
# raw TPCs for 111, 211 and 212 are 0, 1/9 and 2/9. All are safe.
# The original nearest-target winner is 212, which is increasing.
# The fallback must reuse those fits and select 211, closest among {111,211}.
fallback_paths <- rbind(matrix(1L,3,3),
                        matrix(rep(c(2L,1L,1L),each=3),3),
                        matrix(rep(c(2L,1L,2L),each=3),3))
fallback_y <- matrix(0L,9,3)
fallback_y[4,1] <- fallback_y[7,3] <- fallback_y[8,2] <- 1L
fallback_fit <- fixture
fallback_fit$observed$dose <- as.vector(fallback_paths)
fallback_fit$observed$tox <- as.vector(fallback_y)
fallback_fit$counts <- dose_cycle_counts(fallback_fit$observed,3,3)
without <- select_mts(fallback_fit, weights="equal", nearest_tie="higher",
                       screen_min_n=0, no_mts_fallback="none")
with_fallback <- select_mts(fallback_fit, weights="equal", nearest_tie="higher",
                             screen_min_n=0)
check("The no-MTS fixture has only the increasing primary winner 212",
      is.na(without$recommendation) && all(without$per_order$selected=="2,1,2"))
check("Restricted revoting selects safe observed 211, closest within the nonincreasing set",
      with_fallback$fallback_used && with_fallback$recommendation=="2,1,1" &&
      all(with_fallback$per_order$fallback_selected=="2,1,1"))
check("The fallback preserves the original PAVA fits and primary winners",
      identical(without$fitted,with_fallback$fitted) &&
      identical(without$per_order$selected,with_fallback$per_order$selected))
stopped_fixture <- fallback_fit
stopped_fixture$stopped <- TRUE
stopped_mts <- select_mts(stopped_fixture,weights="equal",nearest_tie="higher")
check("The fallback never overrides early stopping",
      is.na(stopped_mts$recommendation) && !stopped_mts$fallback_used)
unsafe_fixture <- fallback_fit
unsafe_fixture$snapshots[[1]]$excluded[,] <- TRUE
unsafe_mts <- select_mts(unsafe_fixture,weights="equal",nearest_tie="higher")
check("The fallback never rescues excluded sequences",
      is.na(unsafe_mts$recommendation) && !unsafe_mts$fallback_used)
increasing_fixture <- fallback_fit
increasing_fixture$observed$dose <- rep(c(2L,1L,2L),each=9)
increasing_fixture$observed$tox <- 0L
increasing_fixture$counts <- dose_cycle_counts(increasing_fixture$observed,3,3)
increasing_mts <- select_mts(increasing_fixture,weights="equal",nearest_tie="higher")
check("No MTS remains when no observed safe nonincreasing sequence exists",
      is.na(increasing_mts$recommendation) && !increasing_mts$fallback_used)

dir.create(lyu_path("results"), showWarnings = FALSE)
write.csv(checks, lyu_path("results", "validation_checks.csv"), row.names = FALSE)
write.csv(fit$assignments, lyu_path("results", "assignments.csv"), row.names = FALSE)
write.csv(fit$observed, lyu_path("results", "observed_outcomes.csv"), row.names = FALSE)
write.csv(fit$continuation, lyu_path("results", "continuation.csv"), row.names = FALSE)
write.csv(totals, lyu_path("results", "cumulative_counts_comparison.csv"), row.names = FALSE)
write.csv(appendix_table, lyu_path("results", "appendix_decision_table.csv"), row.names = FALSE)
write.csv(lower$per_order, lyu_path("results", "mts_ordering_votes.csv"), row.names = FALSE)
print(checks, row.names = FALSE)
cat("\n", nrow(checks), "checks passed.\n", sep = "")
