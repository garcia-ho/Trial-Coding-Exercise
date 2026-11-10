# Table 3 simulation. Run on Unity from the shared project root:
# Rscript --vanilla Lyu_2018/simulate.R 10000 32 2018
# Arguments: trials PER SCENARIO, allocated CPU cores, seed. All six scenarios.
# Only data/table3.rds is written. source() defines functions without simulating.
lyu_dir <- if (file.exists("Lyu_2018/R/mtpi2.R")) "Lyu_2018" else "."
for (name in c("mtpi2", "basyc", "mts")) {
  source(file.path(lyu_dir, "R", paste0(name, ".R")))
}
lyu_scenarios <- function() {
  baseline <- rbind(c(.15, .30, .45), c(.15, .27, .45), c(.08, .16, .24),
                    c(.48, .66, .84), c(.08, .16, .58), c(.15, .46, .53))
  # Scenario 4 gives .84 * 1.1^2 = 1.0164 at D3/C3. The paper's stated
  # multiplier needs a probability cap; record this extra convention explicitly.
  lapply(seq_len(nrow(baseline)), function(s) pmin(outer(baseline[s, ], 1.1^(0:2)), 1))
}

# Enumerate nonincreasing paths and apply the simulation truth definition.
# Safety here uses TRUE p < .35, separately from posterior trial safety.
lyu_true_mts <- function(p, target = .30, ei = c(.25, .35)) {
  paths <- as.matrix(expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:3))
  paths <- paths[apply(paths, 1, function(d) all(diff(d) <= 0)), , drop = FALSE]
  safe <- apply(paths, 1, function(d) all(p[cbind(d, 1:3)] < ei[2]))
  paths <- paths[safe, , drop = FALSE]
  if (!nrow(paths)) return("None")
  q <- apply(paths, 1, function(d) mean(p[cbind(d, 1:3)]))
  eligible <- q > ei[1] & q < ei[2]
  if (!any(eligible)) eligible <- q < target
  if (!any(eligible)) return("None")
  paths <- paths[eligible, , drop = FALSE]
  q <- q[eligible]
  ord <- order(abs(q - target), -rowSums(paths), -paths[, 1], -paths[, 2])
  paste(paths[ord[1], ], collapse = ",")
}

# Independent patient/cycle Bernoulli outcomes are an explicit simulation
# assumption; the paper does not provide a complete longitudinal generator.
simulate_table3 <- function(nsim = 10000L, cores = 32L, seed = 2018L,
                            weights = "equal", nearest_tie = "higher",
                            safety_min_n = 4L, scenarios = 1:6,
                            mts_screen_min_n = 0L, no_mts_fallback = "restricted_revote",
                            global_monitor = "after_first_sequence") {
  start_time <- proc.time()[["elapsed"]]
  stopifnot(length(nsim) == 1, nsim >= 1, nsim <= 100000, nsim %% 1 == 0,
            length(seed) == 1, seed >= 1, seed <= .Machine$integer.max - 700000,
            seed %% 1 == 0, safety_min_n %in% c(0, 3, 4), mts_screen_min_n %in% c(0, 3),
            length(cores) == 1, cores >= 1, cores %% 1 == 0,
            length(scenarios) > 0, all(scenarios %in% 1:6), !anyDuplicated(scenarios))
  # Never use more cores than Slurm allocated to this R process.
  allocated <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = cores))
  cores <- min(as.integer(cores), allocated, length(scenarios) * nsim)
  if (.Platform$OS.type == "windows") cores <- 1L
  weights <- match.arg(weights, c("equal", "patient_cycles"))
  # Restore original PAVA position ties. Keep the existing fallback.
  nearest_tie <- match.arg(nearest_tie, c("lower", "higher"))
  no_mts_fallback <- match.arg(no_mts_fallback, c("restricted_revote", "none"))
  global_monitor <- match.arg(global_monitor,c("continuous","after_first_sequence"))
  truth <- lyu_scenarios()
  target_mts <- vapply(truth, lyu_true_mts, character(1))
  # Strict intratrial n > 3 is implemented as minimum n=4 in the shared
  # safety checks (global stop, DCR and permanent dose exclusion).
  # The final MTS posterior screen retains the appendix's unqualified rule
  # (minimum n=0). Keep these choices separate rather than changing both.
  settings <- basyc_settings(safety_min_n = safety_min_n)
  # Split trial indices into batches, allowing more workers than scenarios.
  # Fixed trial seeds make worker count/scheduling irrelevant.
  tasks <- expand.grid(trial = seq_len(nsim), scenario = scenarios)
  batches <- split(seq_len(nrow(tasks)),
                   rep(seq_len(min(nrow(tasks), 2L * cores)), length.out = nrow(tasks)))
  cat("Running", nrow(tasks), "trials across", cores, "allocated cores.\n")
  simulate_batch <- function(ids) {
    trials <- cells <- paths <- vector("list", length(ids))
    for (r in seq_along(ids)) {
      s <- tasks$scenario[ids[r]]
      b <- tasks$trial[ids[r]]
      # Separate reproducible seeds: extending nsim preserves earlier trials.
      trial_seed <- as.integer(seed + (s - 1L) * 100000L + b)
      set.seed(trial_seed)
      uniforms <- matrix(runif(15L * 3L), nrow = 15L)
      fit <- run_basyc(uniforms, settings = settings, toxicity_prob = truth[[s]],
                        global_monitor=global_monitor)
      mts <- select_mts(fit, weights = weights, nearest_tie = nearest_tie,
                         screen_min_n = mts_screen_min_n,
                         no_mts_fallback = no_mts_fallback)
      selected <- if (is.na(mts$recommendation)) "None" else mts$recommendation
      o <- fit$observed
      trials[[r]] <- data.frame(
        scenario = s, trial = b, seed = trial_seed, selected_mts = selected,
        true_mts = target_mts[s], correct_selection = selected == target_mts[s],
        stopped_early = fit$stopped, fallback_used = mts$fallback_used,
        n_enrolled = length(unique(fit$assignments$patient)),
        n_complete = sum(tabulate(o$patient, nbins = 15L) == 3L),
        patient_cycles = nrow(o), dlts = sum(o$tox), end_week = fit$end_week,
        global_start_week=fit$global_start_week,stop_cycle=fit$stop_event$cycle,
        stop_x=fit$stop_event$x,stop_n=fit$stop_event$n,stop_risk=fit$stop_event$risk)
      z <- expand.grid(dose = 1:3, cycle = 1:3)
      z$scenario <- s
      z$trial <- b
      z$n <- as.vector(fit$counts$n)
      z$dlts <- as.vector(fit$counts$x)
      cells[[r]] <- z
      complete <- which(tabulate(o$patient, nbins = 15L) == 3L)
      paths[[r]] <- data.frame(scenario = integer(), trial = integer(),
                               sequence = character(), dlts = integer())
      if (length(complete)) {
        paths[[r]] <- do.call(rbind, lapply(complete, function(id) {
          patient <- o[o$patient == id, ]
          patient <- patient[order(patient$cycle), ]
          data.frame(scenario = s, trial = b,
                      sequence = paste(patient$dose, collapse = ","), dlts = sum(patient$tox))
        }))
      }
    }
    cat("Finished a batch of", length(ids), "trials.\n")
    list(trials = do.call(rbind, trials), cells = do.call(rbind, cells),
         paths = do.call(rbind, paths))
  }
  batch_results <- parallel::mclapply(batches, simulate_batch, mc.cores = cores,
                                          mc.set.seed = FALSE)
  if (any(vapply(batch_results, inherits, logical(1), "try-error"))) {
    stop("A simulation worker failed; no result file was written.")
  }
  combine <- function(name) {
    z <- do.call(rbind, lapply(batch_results, `[[`, name))
    z <- z[order(z$scenario, z$trial, method = "radix"), , drop = FALSE]
    rownames(z) <- NULL
    z
  }
  truth_table <- do.call(rbind, lapply(scenarios, function(s) {
    z <- expand.grid(dose = 1:3, cycle = 1:3)
    z$scenario <- s
    z$p <- as.vector(truth[[s]])
    z$true_mts <- target_mts[s]
    z
  }))
  run_settings <- data.frame(
    setting = c("nsim_per_scenario", "scenarios", "cores", "seed", "maximum_patients", "cohort_size",
                "cycles", "target", "EI", "xi", "safety_min_n", "mts_screen_min_n", "pava_weights",
                "nearest_tie", "clinical_ties", "no_mts_fallback", "global_monitor",
                "global_exclusion_start", "global_evaluation", "outcome_generator", "probability_cap", "R_version", "created_utc", "host"),
    value = c(nsim, paste(scenarios, collapse = ","), cores, seed, 15, 3, 3, .30, "(0.25,0.35)", .95, safety_min_n, mts_screen_min_n,
              weights, nearest_tie, FALSE, no_mts_fallback,global_monitor,global_monitor,
              if(global_monitor=="continuous")"all assessments" else "before further treatment",
              "independent patient/cycle Bernoulli", "pmin(p,1)",
              R.version.string, format(Sys.time(), tz = "UTC", usetz = TRUE),
              Sys.info()[["nodename"]]))
  result <- list(trials = combine("trials"), cells = combine("cells"),
       paths = combine("paths"), truth = truth_table, settings = run_settings,
       session = capture.output(sessionInfo()))
  result$runtime_seconds <- unname(proc.time()[["elapsed"]] - start_time)
  cat("Total simulation runtime:", round(result$runtime_seconds, 2), "seconds.\n")
  result
}

# Summarize SAVED trials; this does not simulate or call the decision engine.
summarize_table3 <- function(result) {
  nsim <- as.integer(result$settings$value[result$settings$setting == "nsim_per_scenario"])
  labels <- apply(as.matrix(expand.grid(1:3, 1:3, 1:3)), 1, paste, collapse = ",")
  selection <- do.call(rbind, lapply(sort(unique(result$trials$scenario)), function(s) {
    z <- subset(result$trials, scenario == s)
    frequency <- as.integer(table(factor(z$selected_mts, levels = c(labels, "None"))))
    p <- frequency / nsim
    patients <- subset(result$paths, scenario == s)
    data.frame(scenario = s, sequence = c(labels, "None"),
      selection = p, mcse = sqrt(p * (1-p) / nsim),
      mean_dlts = vapply(c(labels, "None"), function(path) {
        y <- patients$dlts[patients$sequence == path]
        if (length(y)) mean(y) else NA_real_
      }, numeric(1)),
      sd_dlts = vapply(c(labels, "None"), function(path) {
        sd(patients$dlts[patients$sequence == path])
      }, numeric(1)),
      allocation = vapply(c(labels, "None"), function(path) {
        sum(patients$sequence == path) / (15 * nsim)
      }, numeric(1)))
  }))
  cells <- aggregate(cbind(n, dlts) ~ scenario + dose + cycle, result$cells, mean)
  early_stop <- aggregate(stopped_early ~ scenario, result$trials, mean)
  list(selection = selection, cells = cells, early_stop = early_stop)
}

# Printed Table 3 rows, transcribed for direct comparisons in the Rmd.
# Selection and allocation entries are proportions, despite the printed (%) label.
table3_reference <- function() {
  selection <- read.table(text = '
scenario sequence selection mean_dlts sd_dlts allocation
1 1,1,1 .242 .544 .671 .304
1 2,1,1 .197 .902 .769 .164
1 2,2,1 .190 .645 .625 .089
1 2,2,2 .131 .363 .481 .034
1 3,2,1 .070 1.173 .800 .063
1 3,1,1 .056 .869 .756 .035
2 1,1,1 .209 .534 .666 .288
2 2,2,1 .185 .629 .627 .090
2 2,1,1 .174 .869 .771 .145
2 2,2,2 .148 .328 .469 .042
2 3,2,1 .082 1.158 .818 .070
2 3,2,2 .065 .874 .682 .032
2 3,1,1 .057 .883 .761 .039
3 3,3,2 .209 .597 .630 .080
3 3,3,3 .204 .289 .453 .044
3 3,2,2 .163 .579 .627 .083
3 2,2,2 .140 .196 .397 .094
3 3,2,1 .066 .975 .782 .068
3 1,1,1 .061 .276 .503 .215
3 2,2,1 .059 .534 .579 .074
4 1,1,1 .008 1.390 .855 .202
4 2,1,1 NA 1.687 .861 .021
5 2,2,2 .397 .195 .396 .129
5 3,2,2 .183 .889 .616 .107
5 2,2,1 .136 .552 .585 .092
5 1,1,1 .080 .276 .502 .214
5 3,2,1 .073 1.250 .765 .071
5 2,1,1 .052 .758 .728 .052
6 1,1,1 .450 .543 .669 .373
6 2,1,1 .283 .987 .743 .225
6 2,2,1 .136 .758 .625 .062', header = TRUE, stringsAsFactors = FALSE)
  # Rows within each scenario are D1 C1-C3, D2 C1-C3, D3 C1-C3.
  counts <- list(
    c(1,6.591,1.638,9.854,1.983,10.844,1.784,5.847,1.486,4.518,1.315,3.649,1.061,2.373,.177,.348,.078,.141),
    c(.936,6.246,1.546,9.369,1.883,10.434,1.543,5.763,1.439,4.818,1.289,3.930,1.247,2.766,.240,.483,.114,.207),
    c(.340,4.325,.518,5.961,.670,6.899,.809,5.075,1.136,6.468,1.210,6.310,1.334,5.572,.664,2.531,.501,1.738),
    c(3.795,7.849,3.118,5.876,1.962,3.400,.527,.800,.039,.054,.012,.016,.018,.021,0,0,0,0),
    c(.354,4.364,.540,6.222,.689,7.133,1.060,6.583,1.313,7.511,1.363,6.977,2.334,4.023,.785,1.223,.584,.836),
    c(1.203,7.997,1.884,11.432,2.205,12.120,2.665,5.810,1.601,3.196,1.351,2.435,.514,.970,.024,.041,.005,.009))
  cells <- do.call(rbind, lapply(1:6, function(s) {
    x <- matrix(counts[[s]], ncol = 2, byrow = TRUE)
    data.frame(scenario = s, dose = rep(1:3, each = 3), cycle = rep(1:3, 3),
                dlts = x[,1], n = x[,2])
  }))
  list(selection = selection, cells = cells,
       early_stop = data.frame(scenario = 1:6,
                              stopped_early = c(3.32,3.75,.43,99.21,.44,4.01)/100))
}

# A source() call from the Rmd stops here. Only Rscript executes this block.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  nsim <- if (length(args) >= 1) as.integer(args[1]) else 10000L
  cores <- if (length(args) >= 2) as.integer(args[2]) else 32L
  seed <- if (length(args) >= 3) as.integer(args[3]) else 2018L
  result <- simulate_table3(nsim = nsim, cores = cores, seed = seed)
  # Replace the one saved dataset only after the whole simulation succeeds.
  output <- file.path(lyu_dir, "data", "table3.rds")
  temporary <- paste0(output, ".tmp")
  saveRDS(result, temporary)
  if (!file.rename(temporary, output)) stop("Could not save ", output)
  cat("Saved", nrow(result$trials), "trials to", output, "\n")
}
