# Functions used by Drug_combo_Boin.Rmd.
# Educational implementation only; not validated for clinical trial conduct.

# Calculate the optimal BOIN escalation and de-escalation boundaries.
# The conventional defaults compare the target phi with phi1 = 0.6 * phi
# (subtherapeutic) and phi2 = 1.4 * phi (overly toxic).
boin_boundaries <- function(target, phi1 = 0.6 * target,
                            phi2 = 1.4 * target) {
  stopifnot(
    length(target) == 1L, target > 0, target < 1,
    phi1 > 0, phi1 < target,
    phi2 > target, phi2 < 1
  )
  lambda_e <- log((1 - phi1) / (1 - target)) /
    log(target * (1 - phi1) / (phi1 * (1 - target)))
  lambda_d <- log((1 - target) / (1 - phi2)) /
    log(phi2 * (1 - target) / (target * (1 - phi2)))
  c(escalate = lambda_e, deescalate = lambda_d)
}

# Return posterior summaries for one dose combination. Each cell has the
# noninformative Jeffreys prior Beta(0.5, 0.5), as in combination BOIN.
combo_posterior_summary <- function(y, n, boundaries, target,
                                    prior_a = 0.5, prior_b = 0.5) {
  shape1 <- prior_a + y
  shape2 <- prior_b + n - y
  c(
    mean = shape1 / (shape1 + shape2),
    desirable = unname(
      pbeta(boundaries["deescalate"], shape1, shape2) -
        pbeta(boundaries["escalate"], shape1, shape2)
    ),
    overdose = unname(1 - pbeta(target, shape1, shape2))
  )
}

# Check that toxicity is nondecreasing in each drug when the other is held
# fixed. This coordinate-wise monotonicity is the partial-order assumption.
is_monotone_combo_surface <- function(probability, tolerance = 1e-12) {
  row_ok <- all(apply(probability, 1L, function(x) {
    all(diff(x) >= -tolerance)
  }))
  column_ok <- all(apply(probability, 2L, function(x) {
    all(diff(x) >= -tolerance)
  }))
  row_ok && column_ok
}

# Validate a true DLT-probability surface before it is used for simulation.
validate_combo_surface <- function(probability) {
  if (!is.matrix(probability) || any(!is.finite(probability)) ||
      any(probability <= 0) || any(probability >= 1)) {
    stop("The true toxicity surface must be a matrix with values in (0, 1).")
  }
  if (!is_monotone_combo_surface(probability)) {
    stop("The true toxicity surface must be nondecreasing in each drug.")
  }
  invisible(TRUE)
}

# Four 4 x 4 scenarios used in the notebook. Interaction is interpreted
# qualitatively relative to an additive reference surface. Directly specifying
# the full joint surface avoids requiring known monotherapy toxicity curves.
make_combo_scenarios <- function() {
  scenarios <- list(
    "Global synergy" = matrix(c(
      0.04, 0.08, 0.15, 0.24,
      0.08, 0.15, 0.26, 0.39,
      0.15, 0.26, 0.40, 0.55,
      0.24, 0.39, 0.55, 0.70
    ), nrow = 4L, byrow = TRUE),
    "Global antagonism" = matrix(c(
      0.04, 0.07, 0.11, 0.16,
      0.07, 0.11, 0.16, 0.22,
      0.11, 0.16, 0.22, 0.28,
      0.16, 0.22, 0.28, 0.35
    ), nrow = 4L, byrow = TRUE),
    "Low synergy / high antagonism" = matrix(c(
      0.06, 0.12, 0.18, 0.23,
      0.12, 0.19, 0.24, 0.29,
      0.18, 0.24, 0.30, 0.36,
      0.23, 0.29, 0.36, 0.44
    ), nrow = 4L, byrow = TRUE),
    "Low antagonism / high synergy" = matrix(c(
      0.02, 0.04, 0.08, 0.15,
      0.04, 0.08, 0.15, 0.27,
      0.08, 0.15, 0.28, 0.45,
      0.15, 0.27, 0.45, 0.68
    ), nrow = 4L, byrow = TRUE)
  )
  lapply(scenarios, function(x) {
    validate_combo_surface(x)
    dimnames(x) <- list(
      paste0("A", seq_len(nrow(x))),
      paste0("B", seq_len(ncol(x)))
    )
    x
  })
}

# Reference surface used only to explain the four interaction patterns. It is
# not used by the dose-assignment algorithm or by the performance metrics.
make_additive_reference <- function() {
  matrix(c(
    0.04, 0.07, 0.12, 0.20,
    0.07, 0.12, 0.20, 0.30,
    0.12, 0.20, 0.30, 0.43,
    0.20, 0.30, 0.43, 0.58
  ), nrow = 4L, byrow = TRUE,
  dimnames = list(paste0("A", 1:4), paste0("B", 1:4)))
}

# Identify all cells whose true DLT probability is closest to target. A
# combination trial can have several true MTDs along an MTD contour.
true_combo_mtds <- function(true_probability, target) {
  distance <- abs(true_probability - target)
  which(abs(distance - min(distance)) < 1e-12, arr.ind = TRUE)
}

# Build an experimental candidate set of up to 2, 3, or 4 combinations.
#
# With candidate_size = 2 this reproduces the standard BOIN neighborhood:
# change one drug by one level. Size 3 adds the northeast (or southwest)
# diagonal. Size 4 adds an anti-diagonal move that changes B in the intended
# direction while changing A in the opposite direction. For example, from
# A2B2 the escalation set is A3B2, A3B3, A2B3, A1B3. This explores a cell whose
# order relative to the current dose is unknown under partial monotonicity.
combo_candidate_set <- function(current, n_a, n_b, direction, eliminated,
                                candidate_size, n) {
  stopifnot(direction %in% c("escalate", "deescalate"),
            candidate_size %in% 2:4)
  offsets <- if (direction == "escalate") {
    rbind(
      c(1L, 0L),   # increase A
      c(1L, 1L),   # increase both
      c(0L, 1L),   # increase B
      c(-1L, 1L)   # anti-diagonal: lower A, increase B
    )
  } else {
    rbind(
      c(-1L, 0L),  # decrease A
      c(-1L, -1L), # decrease both
      c(0L, -1L),  # decrease B
      c(1L, -1L)   # anti-diagonal: increase A, lower B
    )
  }
  keep_rows <- if (candidate_size == 2L) {
    c(1L, 3L)
  } else if (candidate_size == 3L) {
    1:3
  } else {
    1:4
  }
  candidates <- sweep(offsets[keep_rows, , drop = FALSE], 2L, current, "+")
  valid <- candidates[, 1L] >= 1L & candidates[, 1L] <= n_a &
    candidates[, 2L] >= 1L & candidates[, 2L] <= n_b
  candidates <- candidates[valid, , drop = FALSE]
  if (!nrow(candidates)) return(matrix(integer(), ncol = 2L))
  allowed <- !eliminated[candidates]
  unname(candidates[allowed, , drop = FALSE])
}

# Choose the admissible combination with the largest posterior probability of
# lying in the BOIN stay interval. Random selection resolves exact ties.
choose_combo_candidate <- function(candidates, y, n, boundaries, target) {
  if (!nrow(candidates)) return(NULL)
  scores <- apply(candidates, 1L, function(cell) {
    combo_posterior_summary(
      y[cell[1L], cell[2L]], n[cell[1L], cell[2L]], boundaries, target
    )["desirable"]
  })
  best <- which(abs(scores - max(scores)) < 1e-12)
  unname(candidates[sample(best, 1L), ])
}

# Decide whether the high-toxicity observation specified in reliability metric
# C3 has occurred at the current cell. The second clause extends the stated
# >3/6 rule to later cumulative looks as an observed rate greater than 1/2.
irrational_deescalation_trigger <- function(y, n) {
  (n == 3L && y >= 2L) || (n >= 6L && y / n > 0.5)
}

# Simulate one two-drug combination BOIN trial.
run_combo_boin_trial <- function(
    true_probability,
    target = 0.25,
    n_patients = 36L,
    cohort_size = 3L,
    candidate_size = 2L,
    start = c(1L, 1L),
    elimination_cutoff = 0.95,
    seed = NULL) {
  validate_combo_surface(true_probability)
  stopifnot(length(start) == 2L, candidate_size %in% 2:4,
            n_patients > 0, cohort_size > 0)
  if (!is.null(seed)) set.seed(seed)

  n_a <- nrow(true_probability)
  n_b <- ncol(true_probability)
  boundaries <- boin_boundaries(target)
  n <- matrix(0L, n_a, n_b)
  y <- matrix(0L, n_a, n_b)
  eliminated <- matrix(FALSE, n_a, n_b)
  current <- as.integer(start)
  history_rows <- vector("list", ceiling(n_patients / cohort_size))
  patient <- 0L
  cohort <- 0L
  stopped_early <- FALSE

  while (patient < n_patients) {
    cohort <- cohort + 1L
    cohort_n <- min(cohort_size, n_patients - patient)
    outcomes <- rbinom(
      cohort_n, 1L, true_probability[current[1L], current[2L]]
    )
    n[current[1L], current[2L]] <-
      n[current[1L], current[2L]] + cohort_n
    y[current[1L], current[2L]] <-
      y[current[1L], current[2L]] + sum(outcomes)
    n_here <- n[current[1L], current[2L]]
    y_here <- y[current[1L], current[2L]]
    observed <- y_here / n_here
    posterior <- combo_posterior_summary(
      y_here, n_here, boundaries, target
    )

    # Under coordinate-wise monotonicity, eliminating a cell also eliminates
    # all cells northeast of it (at least as high in both agents).
    if (n_here >= 3L && posterior["overdose"] > elimination_cutoff) {
      eliminated[current[1L]:n_a, current[2L]:n_b] <- TRUE
    }

    action <- if (eliminated[current[1L], current[2L]]) {
      "deescalate"
    } else if (observed <= boundaries["escalate"]) {
      "escalate"
    } else if (observed >= boundaries["deescalate"]) {
      "deescalate"
    } else {
      "stay"
    }

    next_combo <- current
    candidate_label <- ""
    if (eliminated[1L, 1L]) {
      action <- "stop: lowest combination too toxic"
      stopped_early <- TRUE
    } else if (action %in% c("escalate", "deescalate")) {
      candidates <- combo_candidate_set(
        current, n_a, n_b, action, eliminated, candidate_size, n
      )
      if (nrow(candidates)) {
        candidate_label <- paste(
          paste0("A", candidates[, 1L], "B", candidates[, 2L]),
          collapse = ", "
        )
      }
      chosen <- choose_combo_candidate(
        candidates, y, n, boundaries, target
      )
      if (is.null(chosen)) {
        if (eliminated[current[1L], current[2L]]) {
          action <- "stop: no admissible lower combination"
          stopped_early <- TRUE
        } else {
          action <- "stay: boundary of dose grid"
        }
      } else {
        next_combo <- chosen
      }
    }

    c3_opportunity <- irrational_deescalation_trigger(y_here, n_here)
    strict_deescalation <-
      all(next_combo <= current) && any(next_combo < current)
    safe_stop <- stopped_early
    c3_failure <- c3_opportunity && !strict_deescalation && !safe_stop

    history_rows[[cohort]] <- data.frame(
      cohort = cohort,
      first_patient = patient + 1L,
      last_patient = patient + cohort_n,
      drug_a = current[1L],
      drug_b = current[2L],
      true_probability = true_probability[current[1L], current[2L]],
      cohort_dlt = sum(outcomes),
      cumulative_n = n_here,
      cumulative_dlt = y_here,
      observed_rate = observed,
      action = action,
      candidates = candidate_label,
      next_a = next_combo[1L],
      next_b = next_combo[2L],
      c3_opportunity = c3_opportunity,
      c3_failure = c3_failure,
      stringsAsFactors = FALSE
    )
    patient <- patient + cohort_n
    if (stopped_early) break
    current <- next_combo
  }

  history <- do.call(rbind, history_rows[seq_len(cohort)])
  posterior_mean <- (y + 0.5) / (n + 1)
  selectable <- n > 0 & !eliminated
  selected <- c(NA_integer_, NA_integer_)
  if (any(selectable)) {
    distance <- abs(posterior_mean - target)
    distance[!selectable] <- Inf
    best <- which(abs(distance - min(distance)) < 1e-12, arr.ind = TRUE)
    selected <- unname(best[sample(seq_len(nrow(best)), 1L), ])
  }

  list(
    history = history,
    n = n,
    y = y,
    posterior_mean = posterior_mean,
    eliminated = eliminated,
    selected = selected,
    boundaries = boundaries,
    target = target,
    design = paste0(candidate_size, " candidates"),
    candidate_size = candidate_size,
    true_probability = true_probability,
    true_mtds = true_combo_mtds(true_probability, target),
    stopped_early = stopped_early
  )
}

# Simulate local parallel exploration. At any escalation for which two standard
# adjacent combinations and 18 patients remain, three cohorts are precommitted
# to each combination (six cohorts total). Outcomes are evaluated after the
# block, then ordinary two-candidate BOIN resumes until another full block is
# possible. This local rule is not the two-stage PC-BOIN design of Lin and Yin.
run_parallel_combo_boin_trial <- function(
    true_probability,
    target = 0.25,
    n_patients = 36L,
    cohort_size = 3L,
    parallel_cohorts_per_combo = 3L,
    start = c(1L, 1L),
    elimination_cutoff = 0.95,
    seed = NULL) {
  validate_combo_surface(true_probability)
  stopifnot(length(start) == 2L, n_patients > 0, cohort_size > 0,
            parallel_cohorts_per_combo >= 1L)
  if (!is.null(seed)) set.seed(seed)

  n_a <- nrow(true_probability)
  n_b <- ncol(true_probability)
  boundaries <- boin_boundaries(target)
  n <- matrix(0L, n_a, n_b)
  y <- matrix(0L, n_a, n_b)
  eliminated <- matrix(FALSE, n_a, n_b)
  current <- as.integer(start)
  history_rows <- vector("list", ceiling(n_patients / cohort_size))
  patient <- 0L
  cohort <- 0L
  stopped_early <- FALSE
  parallel_blocks <- 0L

  while (patient < n_patients) {
    cohort <- cohort + 1L
    cohort_n <- min(cohort_size, n_patients - patient)
    outcomes <- rbinom(
      cohort_n, 1L, true_probability[current[1L], current[2L]]
    )
    n[current[1L], current[2L]] <-
      n[current[1L], current[2L]] + cohort_n
    y[current[1L], current[2L]] <-
      y[current[1L], current[2L]] + sum(outcomes)
    n_here <- n[current[1L], current[2L]]
    y_here <- y[current[1L], current[2L]]
    observed <- y_here / n_here
    posterior <- combo_posterior_summary(y_here, n_here, boundaries, target)

    if (n_here >= 3L && posterior["overdose"] > elimination_cutoff) {
      eliminated[current[1L]:n_a, current[2L]:n_b] <- TRUE
    }
    action <- if (eliminated[current[1L], current[2L]]) {
      "deescalate"
    } else if (observed <= boundaries["escalate"]) {
      "escalate"
    } else if (observed >= boundaries["deescalate"]) {
      "deescalate"
    } else {
      "stay"
    }

    next_combo <- current
    candidates <- matrix(integer(), ncol = 2L)
    if (action %in% c("escalate", "deescalate")) {
      candidates <- combo_candidate_set(
        current, n_a, n_b, action, eliminated, 2L, n
      )
    }
    candidate_label <- if (nrow(candidates)) {
      paste(paste0("A", candidates[, 1L], "B", candidates[, 2L]),
            collapse = ", ")
    } else {
      ""
    }

    block_patients <- 2L * parallel_cohorts_per_combo * cohort_size
    launch_parallel <- action == "escalate" &&
      nrow(candidates) == 2L &&
      (n_patients - (patient + cohort_n)) >= block_patients

    if (eliminated[1L, 1L]) {
      action <- "stop: lowest combination too toxic"
      stopped_early <- TRUE
    } else if (launch_parallel) {
      action <- "launch parallel exploration"
    } else if (action %in% c("escalate", "deescalate")) {
      chosen <- choose_combo_candidate(
        candidates, y, n, boundaries, target
      )
      if (is.null(chosen)) {
        if (eliminated[current[1L], current[2L]]) {
          action <- "stop: no admissible lower combination"
          stopped_early <- TRUE
        } else {
          action <- "stay: boundary of dose grid"
        }
      } else {
        next_combo <- chosen
      }
    }

    c3_opportunity <- irrational_deescalation_trigger(y_here, n_here)
    strict_deescalation <-
      all(next_combo <= current) && any(next_combo < current)
    c3_failure <- c3_opportunity && !strict_deescalation && !stopped_early
    history_rows[[cohort]] <- data.frame(
      cohort = cohort,
      first_patient = patient + 1L,
      last_patient = patient + cohort_n,
      drug_a = current[1L], drug_b = current[2L],
      true_probability = true_probability[current[1L], current[2L]],
      cohort_dlt = sum(outcomes), cumulative_n = n_here,
      cumulative_dlt = y_here, observed_rate = observed,
      action = action, candidates = candidate_label,
      next_a = if (launch_parallel) NA_integer_ else next_combo[1L],
      next_b = if (launch_parallel) NA_integer_ else next_combo[2L],
      c3_opportunity = c3_opportunity,
      c3_failure = c3_failure,
      stringsAsFactors = FALSE
    )
    patient <- patient + cohort_n
    if (stopped_early) break

    if (launch_parallel) {
      parallel_blocks <- parallel_blocks + 1L
      launch_combo <- current
      last_row_for_candidate <- integer(2L)

      # Interleave the two combinations by round to represent concurrent
      # exploration, while intentionally deferring adaptation until all six
      # precommitted cohorts have outcomes.
      for (round_index in seq_len(parallel_cohorts_per_combo)) {
        for (candidate_index in 1:2) {
          cell <- candidates[candidate_index, ]
          cohort <- cohort + 1L
          outcomes <- rbinom(
            cohort_size, 1L, true_probability[cell[1L], cell[2L]]
          )
          n[cell[1L], cell[2L]] <- n[cell[1L], cell[2L]] + cohort_size
          y[cell[1L], cell[2L]] <-
            y[cell[1L], cell[2L]] + sum(outcomes)
          n_cell <- n[cell[1L], cell[2L]]
          y_cell <- y[cell[1L], cell[2L]]
          history_rows[[cohort]] <- data.frame(
            cohort = cohort,
            first_patient = patient + 1L,
            last_patient = patient + cohort_size,
            drug_a = cell[1L], drug_b = cell[2L],
            true_probability = true_probability[cell[1L], cell[2L]],
            cohort_dlt = sum(outcomes), cumulative_n = n_cell,
            cumulative_dlt = y_cell, observed_rate = y_cell / n_cell,
            action = "parallel block: decision deferred",
            candidates = candidate_label,
            next_a = NA_integer_, next_b = NA_integer_,
            c3_opportunity = FALSE, c3_failure = FALSE,
            stringsAsFactors = FALSE
          )
          last_row_for_candidate[candidate_index] <- cohort
          patient <- patient + cohort_size
        }
      }

      # Apply overdose elimination only after the precommitted parallel block.
      for (candidate_index in 1:2) {
        cell <- candidates[candidate_index, ]
        post <- combo_posterior_summary(
          y[cell[1L], cell[2L]], n[cell[1L], cell[2L]],
          boundaries, target
        )
        if (post["overdose"] > elimination_cutoff) {
          eliminated[cell[1L]:n_a, cell[2L]:n_b] <- TRUE
        }
      }

      available <- !eliminated[cbind(candidates[, 1L], candidates[, 2L])]
      observed_candidates <- vapply(1:2, function(i) {
        y[candidates[i, 1L], candidates[i, 2L]] /
          n[candidates[i, 1L], candidates[i, 2L]]
      }, numeric(1L))
      acceptable <- available &
        observed_candidates < boundaries["deescalate"]
      if (any(acceptable)) {
        current <- choose_combo_candidate(
          candidates[acceptable, , drop = FALSE],
          y, n, boundaries, target
        )
      } else {
        # Both explored combinations look too toxic (or were eliminated), so
        # return to the lower launch combination.
        current <- launch_combo
      }

      # C3 is assessed at the end of the block, the first point at which this
      # parallel strategy permits an adaptive response.
      for (candidate_index in 1:2) {
        cell <- candidates[candidate_index, ]
        row_index <- last_row_for_candidate[candidate_index]
        trigger <- irrational_deescalation_trigger(
          y[cell[1L], cell[2L]], n[cell[1L], cell[2L]]
        )
        strict_lower <- all(current <= cell) && any(current < cell)
        history_rows[[row_index]]$c3_opportunity <- trigger
        history_rows[[row_index]]$c3_failure <- trigger && !strict_lower
      }
      history_rows[[max(last_row_for_candidate)]]$action <-
        "parallel block complete; choose next"
      history_rows[[max(last_row_for_candidate)]]$next_a <- current[1L]
      history_rows[[max(last_row_for_candidate)]]$next_b <- current[2L]
      next
    }

    current <- next_combo
  }

  history <- do.call(rbind, history_rows[seq_len(cohort)])
  posterior_mean <- (y + 0.5) / (n + 1)
  selectable <- n > 0 & !eliminated
  selected <- c(NA_integer_, NA_integer_)
  if (any(selectable)) {
    distance <- abs(posterior_mean - target)
    distance[!selectable] <- Inf
    best <- which(abs(distance - min(distance)) < 1e-12, arr.ind = TRUE)
    selected <- unname(best[sample(seq_len(nrow(best)), 1L), ])
  }

  list(
    history = history,
    n = n,
    y = y,
    posterior_mean = posterior_mean,
    eliminated = eliminated,
    selected = selected,
    boundaries = boundaries,
    target = target,
    design = paste0("Local parallel ", parallel_cohorts_per_combo, "+",
                    parallel_cohorts_per_combo),
    candidate_size = NA_integer_,
    true_probability = true_probability,
    true_mtds = true_combo_mtds(true_probability, target),
    stopped_early = stopped_early,
    parallel_blocks = parallel_blocks,
    parallel_used = parallel_blocks > 0L
  )
}

# Make stable column names for the patient allocation at every grid cell.
allocation_column_names <- function(n_a, n_b) {
  as.vector(outer(
    paste0("A", seq_len(n_a)), paste0("B", seq_len(n_b)),
    function(a, b) paste0("alloc_", a, "_", b)
  ))
}

# Run all scenarios for each experimental candidate-set size. The returned
# data contain one row per trial and enough cell-level information to produce
# selection and mean-allocation bubble maps.
simulate_combo_boin_scenarios <- function(
    scenarios,
    candidate_sizes = 2:4,
    include_parallel = TRUE,
    parallel_cohorts_per_combo = 3L,
    target = 0.25,
    toxic_threshold = 0.33,
    n_sim = 500L,
    n_patients = 36L,
    cohort_size = 3L,
    seed = 20261017L) {
  stopifnot(length(scenarios) > 0L, all(candidate_sizes %in% 2:4),
            toxic_threshold > target, n_sim >= 1L)
  invisible(lapply(scenarios, validate_combo_surface))
  set.seed(seed)
  design_specs <- as.list(candidate_sizes)
  if (include_parallel) design_specs <- c(design_specs, list("parallel"))
  n_records <- length(scenarios) * length(design_specs) * n_sim
  records <- vector("list", n_records)
  index <- 0L

  for (scenario_name in names(scenarios)) {
    truth <- scenarios[[scenario_name]]
    true_mtds <- true_combo_mtds(truth, target)
    allocation_names <- allocation_column_names(nrow(truth), ncol(truth))
    for (design_spec in design_specs) {
      for (simulation in seq_len(n_sim)) {
        trial <- if (identical(design_spec, "parallel")) {
          run_parallel_combo_boin_trial(
            true_probability = truth,
            target = target,
            n_patients = n_patients,
            cohort_size = cohort_size,
            parallel_cohorts_per_combo = parallel_cohorts_per_combo
          )
        } else {
          run_combo_boin_trial(
            true_probability = truth,
            target = target,
            n_patients = n_patients,
            cohort_size = cohort_size,
            candidate_size = design_spec
          )
        }
        h <- trial$history
        selected <- trial$selected
        selected_probability <- if (anyNA(selected)) {
          NA_real_
        } else {
          truth[selected[1L], selected[2L]]
        }
        correct <- !anyNA(selected) && any(
          true_mtds[, 1L] == selected[1L] &
            true_mtds[, 2L] == selected[2L]
        )
        actual_n <- sum(trial$n)
        mtd_patient_count <- sum(trial$n[true_mtds])
        toxic_patient_count <- sum(trial$n[truth > toxic_threshold])
        above_mtd_patient_count <- sum(trial$n[truth > target])
        allocations <- as.numeric(trial$n) / actual_n
        names(allocations) <- allocation_names

        index <- index + 1L
        records[[index]] <- cbind(data.frame(
          scenario = scenario_name,
          design = trial$design,
          candidate_size = trial$candidate_size,
          simulation = simulation,
          selected_a = selected[1L],
          selected_b = selected[2L],
          selected_true_probability = selected_probability,
          correct_selection = correct,
          b1_toxic_selection = !is.na(selected_probability) &&
            selected_probability > toxic_threshold,
          b2_toxic_allocation = toxic_patient_count / actual_n,
          c1_overdose_risk = above_mtd_patient_count / actual_n > 0.50,
          c2_poor_allocation = mtd_patient_count < 6L,
          c3_opportunities = sum(h$c3_opportunity),
          c3_failures = sum(h$c3_failure),
          dlt_rate = sum(trial$y) / actual_n,
          mtd_allocation = mtd_patient_count / actual_n,
          actual_n = actual_n,
          stopped_early = trial$stopped_early,
          stringsAsFactors = FALSE
        ), as.data.frame(as.list(allocations), check.names = FALSE))
      }
    }
  }
  do.call(rbind, records)
}

# Summarize the requested operating characteristics. C3 is event-based: its
# denominator is the total number of qualifying high-toxicity decision points.
summarize_combo_boin <- function(results) {
  group_id <- interaction(
    results$scenario, results$design, drop = TRUE, lex.order = TRUE
  )
  groups <- split(results, group_id)
  answer <- do.call(rbind, lapply(groups, function(x) {
    c3_denominator <- sum(x$c3_opportunities)
    data.frame(
      scenario = x$scenario[1L],
      design = x$design[1L],
      candidate_size = x$candidate_size[1L],
      PCS = mean(x$correct_selection),
      B1_toxic_selection = mean(x$b1_toxic_selection),
      B2_toxic_allocation = mean(x$b2_toxic_allocation),
      C1_overdose_risk = mean(x$c1_overdose_risk),
      C2_poor_allocation = mean(x$c2_poor_allocation),
      C3_irrational_assignment = if (c3_denominator > 0L) {
        sum(x$c3_failures) / c3_denominator
      } else {
        NA_real_
      },
      C3_opportunities = c3_denominator,
      mean_MTD_allocation = mean(x$mtd_allocation),
      mean_DLT_rate = mean(x$dlt_rate),
      early_stop_probability = mean(x$stopped_early),
      mean_actual_n = mean(x$actual_n)
    )
  }))
  rownames(answer) <- NULL
  design_order <- c("2 candidates", "3 candidates", "4 candidates",
                    "Local parallel 3+3")
  answer[order(answer$scenario, match(answer$design, design_order)), ]
}

# A restrained, color-blind-friendly palette used consistently across figures.
combo_series_colors <- function() {
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
}

# Show a true two-drug toxicity surface as four dose-response profiles. Each
# line fixes Drug A and moves across Drug B. Large open rings identify cells on
# the true MTD contour; horizontal lines show phi and the toxic threshold.
plot_combo_scenario_lines <- function(probability, target, toxic_threshold,
                                      main) {
  colors <- combo_series_colors()
  shapes <- c(16L, 17L, 15L, 18L)
  n_a <- nrow(probability)
  n_b <- ncol(probability)
  plot(
    NA, xlim = c(0.8, n_b + 0.2), ylim = c(0, max(0.75, probability)),
    xaxt = "n", yaxt = "n", xlab = "Drug B dose level",
    ylab = "True DLT probability",
    main = main, cex.lab = 1.05, cex.main = 1.1
  )
  axis(1L, at = seq_len(n_b), labels = paste0("B", seq_len(n_b)),
       cex.axis = 1.0)
  axis(2L, las = 1L, cex.axis = 0.95)
  abline(h = pretty(c(0, max(0.75, probability))), col = "#E5E7EB",
         lwd = 1)
  abline(h = target, col = "#374151", lwd = 2)
  abline(h = toxic_threshold, col = "#B91C1C", lwd = 2, lty = 2)
  for (a in seq_len(n_a)) {
    lines(seq_len(n_b), probability[a, ], type = "b", lwd = 2.2,
          pch = shapes[a], cex = 1.15, col = colors[a])
  }
  mtd <- true_combo_mtds(probability, target)
  points(
    mtd[, 2L], probability[mtd], pch = 21L, cex = 2.0,
    lwd = 2.2, bg = NA, col = "#111827"
  )
  legend(
    "topleft", legend = paste0("A", seq_len(n_a)),
    col = colors, pch = shapes, lty = 1L, lwd = 2,
    ncol = 2L, bty = "n", cex = 0.82
  )
  legend(
    "bottomright", legend = c("Target 25%", "Toxic >33%", "True MTD"),
    col = c("#374151", "#B91C1C", "#111827"),
    lty = c(1L, 2L, NA), pch = c(NA, NA, 21L),
    pt.bg = c(NA, NA, NA), lwd = c(2, 2, NA),
    bty = "n", cex = 0.76
  )
}

# Bubble map for a 4 x 4 probability or allocation surface. Circle size and
# direct percentage labels carry the value; the pale fill is only supportive,
# so the plot remains readable without relying on a dark color scale.
plot_combo_bubble_surface <- function(probability, main, true_mtds = NULL,
                                      scale_max = 0.50) {
  n_a <- nrow(probability)
  n_b <- ncol(probability)
  plot(
    NA, xlim = c(0.5, n_b + 0.5), ylim = c(0.5, n_a + 0.5),
    xaxt = "n", yaxt = "n", xlab = "Drug B", ylab = "Drug A",
    main = main, asp = 1, cex.lab = 1.0, cex.main = 0.95
  )
  axis(1L, at = seq_len(n_b), labels = paste0("B", seq_len(n_b)),
       cex.axis = 0.95)
  axis(2L, at = seq_len(n_a), labels = paste0("A", seq_len(n_a)),
       cex.axis = 0.95, las = 1L)
  abline(v = seq_len(n_b), h = seq_len(n_a), col = "#E5E7EB", lwd = 1)
  grid <- expand.grid(drug_a = seq_len(n_a), drug_b = seq_len(n_b))
  values <- probability[cbind(grid$drug_a, grid$drug_b)]
  positive <- values >= 0.005
  point_size <- 1.0 + 3.4 * sqrt(pmin(values[positive] / scale_max, 1))
  points(
    grid$drug_b[positive], grid$drug_a[positive], pch = 21L,
    cex = point_size,
    bg = grDevices::adjustcolor("#56B4E9", alpha.f = 0.35),
    col = "#0072B2", lwd = 1.3
  )
  text(
    grid$drug_b[positive], grid$drug_a[positive],
    labels = paste0(round(100 * values[positive]), "%"),
    cex = 0.92, col = "#111827", font = 2L
  )
  if (!is.null(true_mtds)) {
    apply(true_mtds, 1L, function(cell) {
      rect(cell[2L] - 0.43, cell[1L] - 0.43,
           cell[2L] + 0.43, cell[1L] + 0.43,
           border = "#D55E00", lwd = 2.4)
    })
  }
  box(col = "#9CA3AF")
}

# Plot the six operating characteristics as line-and-dot small multiples. A
# shared mapping is used: color and point shape both identify the scenario.
plot_combo_operating_dots <- function(summary_table) {
  metric_names <- c(
    "PCS", "B1_toxic_selection", "B2_toxic_allocation",
    "C1_overdose_risk", "C2_poor_allocation",
    "C3_irrational_assignment"
  )
  metric_labels <- c(
    "Accuracy: PCS", "B1: toxic selection", "B2: toxic allocation",
    "C1: >50% above target", "C2: <6 at true MTD",
    "C3: irrational assignment"
  )
  scenario_levels <- unique(summary_table$scenario)
  design_levels <- unique(summary_table$design)
  colors <- combo_series_colors()
  shapes <- c(16L, 17L, 15L, 18L)
  old <- par(mfrow = c(3, 2), mar = c(5.2, 4.8, 3, 1),
             oma = c(0, 0, 1, 0))
  on.exit(par(old))
  for (k in seq_along(metric_names)) {
    values <- summary_table[[metric_names[k]]]
    upper <- min(1, max(0.10, values, na.rm = TRUE) * 1.18)
    plot(
      NA, xlim = c(0.75, length(design_levels) + 0.25),
      ylim = c(0, upper), xaxt = "n", yaxt = "n", xlab = "Design",
      ylab = "Probability / proportion", main = metric_labels[k],
      cex.lab = 1.1, cex.main = 1.18
    )
    axis(1L, at = seq_along(design_levels),
         labels = c("Size 2", "Size 3", "Size 4", "Local\n3+3"),
         cex.axis = 1.0)
    axis(2L, las = 1L, cex.axis = 1.0)
    abline(h = pretty(c(0, upper)), col = "#E5E7EB", lwd = 1)
    for (i in seq_along(scenario_levels)) {
      row <- summary_table$scenario == scenario_levels[i]
      y_values <- summary_table[[metric_names[k]]][row]
      lines(seq_along(design_levels), y_values, type = "b",
            col = colors[i], pch = shapes[i], lwd = 2, cex = 1.1)
    }
    if (k == 1L) {
      legend(
        "topright", legend = scenario_levels, col = colors,
        pch = shapes, lty = 1L, lwd = 2, bty = "n", cex = 0.82
      )
    }
  }
}

# Plot the four designs for one scenario in a two-column bubble-map page.
plot_combo_selection_bubbles <- function(results, scenarios, target,
                                         scenario_name) {
  design_levels <- unique(results$design)
  old <- par(
    mfrow = c(2, 2), mar = c(4.0, 4.2, 3.2, 1), oma = c(0, 0, 2.2, 0)
  )
  on.exit(par(old))
  truth <- scenarios[[scenario_name]]
  for (design_name in design_levels) {
    group <- results$scenario == scenario_name &
      results$design == design_name
    selected <- results[group & !is.na(results$selected_a), ]
    probability <- matrix(0, nrow(truth), ncol(truth))
    if (nrow(selected)) {
      counts <- table(selected$selected_a, selected$selected_b)
      for (a in rownames(counts)) for (b in colnames(counts)) {
        probability[as.integer(a), as.integer(b)] <-
          counts[a, b] / sum(group)
      }
    }
    plot_combo_bubble_surface(
      probability, main = design_name,
      true_mtds = true_combo_mtds(truth, target), scale_max = 0.50
    )
  }
  mtext(
    paste0(scenario_name,
           ": selection probability; orange boxes mark true MTDs"),
    outer = TRUE, cex = 1.18, font = 2L
  )
}

# Plot mean patient allocation for one scenario in a two-column page.
plot_combo_allocation_bubbles <- function(results, scenarios, scenario_name) {
  design_levels <- unique(results$design)
  allocation_names <- grep("^alloc_A[0-9]+_B[0-9]+$", names(results),
                           value = TRUE)
  old <- par(
    mfrow = c(2, 2), mar = c(4.0, 4.2, 3.2, 1), oma = c(0, 0, 2.2, 0)
  )
  on.exit(par(old))
  truth <- scenarios[[scenario_name]]
  for (design_name in design_levels) {
    group <- results$scenario == scenario_name &
      results$design == design_name
    mean_allocation <- colMeans(results[group, allocation_names,
                                          drop = FALSE])
    allocation_matrix <- matrix(
      mean_allocation, nrow = nrow(truth), ncol = ncol(truth)
    )
    plot_combo_bubble_surface(
      allocation_matrix, main = design_name, scale_max = 0.50
    )
  }
  mtext(paste0(scenario_name, ": mean patient-allocation proportion"),
        outer = TRUE, cex = 1.18, font = 2L)
}

# Plot the sequence of visited combinations for one simulated trial. Point
# color identifies whether at least one DLT occurred in that cohort.
plot_combo_trial_path <- function(trial, main = NULL) {
  h <- trial$history
  if (is.null(main)) {
    main <- trial$design
  }

  # Separate repeated visits to the same combination by a small deterministic
  # offset. This keeps cohort numbers legible without moving a point far enough
  # to suggest that a different dose combination was used.
  cell_id <- paste(h$drug_a, h$drug_b, sep = "-")
  visit_order <- ave(seq_len(nrow(h)), cell_id, FUN = seq_along)
  visit_total <- ave(seq_len(nrow(h)), cell_id, FUN = length)
  angle <- 2 * pi * (visit_order - 1) / pmax(visit_total, 1)
  radius <- ifelse(visit_total > 1, 0.16, 0)
  plot_a <- h$drug_a + radius * sin(angle)
  plot_b <- h$drug_b + radius * cos(angle)

  plot(
    plot_b, plot_a, type = "o", lwd = 2.2, pch = 21L, cex = 1.55,
    bg = ifelse(h$cohort_dlt > 0, "#FDBA74", "#BAE6FD"),
    col = "#374151",
    xlim = c(0.7, ncol(trial$true_probability) + 0.3),
    ylim = c(0.7, nrow(trial$true_probability) + 0.3),
    xlab = "Drug B level", ylab = "Drug A level", main = main,
    cex.lab = 1.15, cex.main = 1.15,
    xaxt = "n", yaxt = "n"
  )
  axis(1L, at = seq_len(ncol(trial$true_probability)),
       labels = paste0("B", seq_len(ncol(trial$true_probability))),
       cex.axis = 1.05)
  axis(2L, at = seq_len(nrow(trial$true_probability)),
       labels = paste0("A", seq_len(nrow(trial$true_probability))),
       cex.axis = 1.05, las = 1L)
  text(plot_b, plot_a, labels = h$cohort, cex = 0.80,
       col = "#111827", font = 2L)
  if (!anyNA(trial$selected)) {
    points(trial$selected[2L], trial$selected[1L],
           pch = 4L, cex = 2, lwd = 2.5, col = "dodgerblue4")
  }
}
