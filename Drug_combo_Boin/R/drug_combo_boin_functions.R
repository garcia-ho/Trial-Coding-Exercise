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

# Enforce the coordinate-wise toxicity ordering at the end of a trial. The
# regression is fitted only to visited cells. It minimizes weighted squared
# error in the Jeffreys posterior means subject to theta(a, b) <= theta(a', b')
# whenever a <= a' and b <= b'. This prevents an isolated DLT at a very low
# combination from making that cell look more toxic than observed higher cells.
combo_isotonic_estimates <- function(y, n, prior_a = 0.5, prior_b = 0.5) {
  stopifnot(is.matrix(y), is.matrix(n), identical(dim(y), dim(n)))
  observed <- which(n > 0L, arr.ind = TRUE)
  estimate <- matrix(NA_real_, nrow(n), ncol(n), dimnames = dimnames(n))
  if (!nrow(observed)) return(estimate)

  values <- (y[observed] + prior_a) /
    (n[observed] + prior_a + prior_b)
  weights <- n[observed] + prior_a + prior_b
  if (nrow(observed) == 1L) {
    estimate[observed] <- values
    return(estimate)
  }

  # Each row of ui represents one inequality ui %*% theta >= ci.
  constraint_rows <- list()
  constraint_index <- 0L
  for (lower in seq_len(nrow(observed) - 1L)) {
    for (higher in (lower + 1L):nrow(observed)) {
      first_below_second <-
        all(observed[lower, ] <= observed[higher, ]) &&
        any(observed[lower, ] < observed[higher, ])
      second_below_first <-
        all(observed[higher, ] <= observed[lower, ]) &&
        any(observed[higher, ] < observed[lower, ])
      if (first_below_second || second_below_first) {
        low <- if (first_below_second) lower else higher
        high <- if (first_below_second) higher else lower
        row <- numeric(nrow(observed))
        row[low] <- -1
        row[high] <- 1
        constraint_index <- constraint_index + 1L
        constraint_rows[[constraint_index]] <- row
      }
    }
  }

  order_constraints <- if (length(constraint_rows)) {
    do.call(rbind, constraint_rows)
  } else {
    matrix(numeric(), nrow = 0L, ncol = nrow(observed))
  }
  # The bounds keep the numerical optimizer in the probability scale.
  ui <- rbind(order_constraints, diag(nrow(observed)),
              -diag(nrow(observed)))
  ci <- c(rep(0, nrow(order_constraints)),
          rep(1e-8, nrow(observed)), rep(-1 + 1e-8, nrow(observed)))

  # This increasing surface is strictly feasible for all order constraints.
  start <- 0.05 + 0.10 *
    ((observed[, 1L] - 1L) + (observed[, 2L] - 1L))
  objective <- function(theta) sum(weights * (theta - values)^2)
  gradient <- function(theta) 2 * weights * (theta - values)
  fitted <- tryCatch(
    constrOptim(start, objective, gradient, ui = ui, ci = ci,
                method = "BFGS", control = list(reltol = 1e-10))$par,
    error = function(e) values
  )
  estimate[observed] <- fitted
  estimate
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

# Return only the two axial neighbors used by local parallel 1+1. This helper
# deliberately excludes diagonal and anti-diagonal moves:
#   escalate from AjBk -> A(j+1)Bk and AjB(k+1)
#   de-escalate from AjBk -> A(j-1)Bk and AjB(k-1)
# Edge cells and combinations already eliminated by the safety rule are removed.
local_parallel_candidate_set <- function(current, n_a, n_b, direction,
                                         eliminated) {
  stopifnot(direction %in% c("escalate", "deescalate"))
  offsets <- if (direction == "escalate") {
    rbind(c(1L, 0L), c(0L, 1L))
  } else {
    rbind(c(-1L, 0L), c(0L, -1L))
  }
  candidates <- sweep(offsets, 2L, current, "+")
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
    n_patients = 48L,
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
  isotonic_mean <- combo_isotonic_estimates(y, n)
  selectable <- n > 0 & !eliminated
  selected <- c(NA_integer_, NA_integer_)
  if (any(selectable)) {
    distance <- abs(isotonic_mean - target)
    distance[!selectable] <- Inf
    best <- which(abs(distance - min(distance)) < 1e-8, arr.ind = TRUE)
    best_n <- n[best]
    best <- best[best_n == max(best_n), , drop = FALSE]
    selected <- unname(best[sample(seq_len(nrow(best)), 1L), ])
  }

  list(
    history = history,
    n = n,
    y = y,
    posterior_mean = posterior_mean,
    isotonic_mean = isotonic_mean,
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

# Simulate local parallel exploration. Whenever a BOIN escalation or
# de-escalation decision has two admissible size-2 neighbors, one cohort is
# sent to each neighbor. The posterior comparison is made after both cohorts
# finish, and the better admissible neighbor becomes the next current dose.
# This repeated two-branch search is an experimental comparator, not a
# published BOIN design.
run_parallel_combo_boin_trial <- function(
    true_probability,
    target = 0.25,
    n_patients = 48L,
    cohort_size = 3L,
    parallel_cohorts_per_combo = 1L,
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
      candidates <- local_parallel_candidate_set(
        current, n_a, n_b, action, eliminated
      )
    }
    candidate_label <- if (nrow(candidates)) {
      paste(paste0("A", candidates[, 1L], "B", candidates[, 2L]),
            collapse = ", ")
    } else {
      ""
    }

    block_patients <- 2L * parallel_cohorts_per_combo * cohort_size
    launch_direction <- action
    launch_parallel <- action %in% c("escalate", "deescalate") &&
      nrow(candidates) == 2L &&
      (n_patients - (patient + cohort_n)) >= block_patients

    if (eliminated[1L, 1L]) {
      action <- "stop: lowest combination too toxic"
      stopped_early <- TRUE
    } else if (launch_parallel) {
      action <- paste("launch parallel", launch_direction, "exploration")
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
    if (launch_parallel && launch_direction == "deescalate") {
      strict_deescalation <- all(apply(candidates, 1L, function(cell) {
        all(cell <= current) && any(cell < current)
      }))
    }
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
      parallel_block = if (launch_parallel) {
        parallel_blocks + 1L
      } else {
        NA_integer_
      },
      parallel_role = if (launch_parallel) "launch" else "sequential",
      parallel_parent_a = if (launch_parallel) current[1L] else NA_integer_,
      parallel_parent_b = if (launch_parallel) current[2L] else NA_integer_,
      stringsAsFactors = FALSE
    )
    patient <- patient + cohort_n
    if (stopped_early) break

    if (launch_parallel) {
      parallel_blocks <- parallel_blocks + 1L
      launch_combo <- current
      last_row_for_candidate <- integer(2L)

      # Interleave the candidates by round to represent concurrent exploration.
      # With the default 1+1 rule, this loop assigns one cohort to each branch
      # and then immediately compares their posteriors.
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
            parallel_block = parallel_blocks,
            parallel_role = paste0("branch ", candidate_index),
            parallel_parent_a = launch_combo[1L],
            parallel_parent_b = launch_combo[2L],
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
      observed_candidates <- vapply(seq_len(nrow(candidates)), function(i) {
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
      } else if (launch_direction == "escalate" &&
                 !eliminated[launch_combo[1L], launch_combo[2L]]) {
        # If both higher branches are eliminated, return to the safe launch
        # combination rather than forcing further escalation.
        current <- launch_combo
      } else {
        # After a de-escalation block, both candidate branches can in principle
        # be eliminated. Fall back to any previously visited, admissible cell
        # that is strictly coordinate-wise lower than the launch combination.
        lower <- which(
          n > 0L & !eliminated &
            row(n) <= launch_combo[1L] & col(n) <= launch_combo[2L] &
            (row(n) < launch_combo[1L] | col(n) < launch_combo[2L]),
          arr.ind = TRUE
        )
        if (nrow(lower)) {
          current <- choose_combo_candidate(
            lower, y, n, boundaries, target
          )
        } else {
          stopped_early <- TRUE
        }
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
        history_rows[[row_index]]$c3_failure <-
          trigger && !strict_lower && !stopped_early
      }
      last_block_row <- max(last_row_for_candidate)
      history_rows[[last_block_row]]$action <- if (stopped_early) {
        "parallel block complete; no admissible lower combination"
      } else {
        "parallel block complete; choose next"
      }
      history_rows[[last_block_row]]$next_a <- if (stopped_early) {
        NA_integer_
      } else {
        current[1L]
      }
      history_rows[[last_block_row]]$next_b <- if (stopped_early) {
        NA_integer_
      } else {
        current[2L]
      }
      if (stopped_early) break
      next
    }

    current <- next_combo
  }

  history <- do.call(rbind, history_rows[seq_len(cohort)])
  posterior_mean <- (y + 0.5) / (n + 1)
  isotonic_mean <- combo_isotonic_estimates(y, n)
  selectable <- n > 0 & !eliminated
  selected <- c(NA_integer_, NA_integer_)
  if (any(selectable)) {
    distance <- abs(isotonic_mean - target)
    distance[!selectable] <- Inf
    best <- which(abs(distance - min(distance)) < 1e-8, arr.ind = TRUE)
    best_n <- n[best]
    best <- best[best_n == max(best_n), , drop = FALSE]
    selected <- unname(best[sample(seq_len(nrow(best)), 1L), ])
  }

  list(
    history = history,
    n = n,
    y = y,
    posterior_mean = posterior_mean,
    isotonic_mean = isotonic_mean,
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

# Run all scenarios for each experimental candidate-set size. Independent
# scenario-design cells can be sent to forked workers on macOS/Linux. Each cell
# has its own deterministic seed, so results do not depend on the worker count.
# The returned data contain one row per trial and enough cell-level information
# to produce selection and mean-allocation bubble maps.
simulate_combo_boin_scenarios <- function(
    scenarios,
    candidate_sizes = 2:4,
    include_parallel = TRUE,
    parallel_cohorts_per_combo = 1L,
    target = 0.25,
    toxic_threshold = 0.33,
    n_sim = 1000L,
    n_patients = 48L,
    cohort_size = 3L,
    seed = 20261017L,
    n_cores = 1L) {
  stopifnot(length(scenarios) > 0L, all(candidate_sizes %in% 2:4),
            toxic_threshold > target, n_sim >= 1L)
  invisible(lapply(scenarios, validate_combo_surface))
  n_cores <- max(1L, as.integer(n_cores))
  design_specs <- as.list(candidate_sizes)
  if (include_parallel) design_specs <- c(design_specs, list("parallel"))
  task_grid <- expand.grid(
    design_index = seq_along(design_specs),
    scenario_index = seq_along(scenarios),
    KEEP.OUT.ATTRS = FALSE
  )

  run_one_cell <- function(task_index) {
    scenario_index <- task_grid$scenario_index[task_index]
    design_index <- task_grid$design_index[task_index]
    scenario_name <- names(scenarios)[scenario_index]
    truth <- scenarios[[scenario_name]]
    true_mtds <- true_combo_mtds(truth, target)
    allocation_names <- allocation_column_names(nrow(truth), ncol(truth))
    design_spec <- design_specs[[design_index]]
    # A stable scenario-design seed lets one method be rerun without changing
    # the random streams for the methods retained from a previous simulation.
    design_seed_index <- if (identical(design_spec, "parallel")) {
      4L
    } else {
      match(as.integer(design_spec), 2:4)
    }
    cell_seed_index <- (scenario_index - 1L) * 4L + design_seed_index
    set.seed(as.integer(seed + 100003L * cell_seed_index))
    records <- vector("list", n_sim)

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

      records[[simulation]] <- cbind(data.frame(
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
    do.call(rbind, records)
  }

  task_indices <- seq_len(nrow(task_grid))
  use_forking <- .Platform$OS.type != "windows" && n_cores > 1L
  pieces <- if (use_forking) {
    parallel::mclapply(
      task_indices, run_one_cell,
      mc.cores = min(n_cores, length(task_indices)),
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(task_indices, run_one_cell)
  }
  answer <- do.call(rbind, pieces)
  rownames(answer) <- NULL
  answer
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
  design_rank <- match(
    answer$design, c("2 candidates", "3 candidates", "4 candidates")
  )
  design_rank[is.na(design_rank)] <- 4L
  answer[order(answer$scenario, design_rank), ]
}

# Repeat the full scenario-by-design experiment over several maximum sample
# sizes. On macOS/Linux, each sample size is assigned to a forked base-R
# worker. A separate, deterministic seed is used for every sample size, so the
# results do not depend on how many workers are requested. Windows uses the
# same code sequentially because base R does not support forked workers there.
simulate_combo_boin_sample_size_grid <- function(
    scenarios,
    sample_sizes = seq(48L, 90L, by = 6L),
    candidate_sizes = 2:4,
    include_parallel = TRUE,
    parallel_cohorts_per_combo = 1L,
    target = 0.25,
    toxic_threshold = 0.33,
    n_sim = 1000L,
    cohort_size = 3L,
    seed = 20261117L,
    n_cores = 6L) {
  sample_sizes <- as.integer(sample_sizes)
  n_cores <- max(1L, as.integer(n_cores))
  stopifnot(
    length(sample_sizes) > 0L,
    all(sample_sizes > 0L),
    all(sample_sizes %% cohort_size == 0L),
    n_sim >= 1L
  )

  run_one_size <- function(i) {
    trial_results <- simulate_combo_boin_scenarios(
      scenarios = scenarios,
      candidate_sizes = candidate_sizes,
      include_parallel = include_parallel,
      parallel_cohorts_per_combo = parallel_cohorts_per_combo,
      target = target,
      toxic_threshold = toxic_threshold,
      n_sim = n_sim,
      n_patients = sample_sizes[i],
      cohort_size = cohort_size,
      seed = as.integer(seed + 100003L * i)
    )
    answer <- summarize_combo_boin(trial_results)
    answer$maximum_sample_size <- sample_sizes[i]
    answer$replicates_per_cell <- n_sim
    answer
  }

  indices <- seq_along(sample_sizes)
  use_forking <- .Platform$OS.type != "windows" && n_cores > 1L
  pieces <- if (use_forking) {
    parallel::mclapply(
      indices, run_one_size,
      mc.cores = min(n_cores, length(indices)),
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(indices, run_one_size)
  }

  answer <- do.call(rbind, pieces)
  rownames(answer) <- NULL
  design_rank <- match(
    answer$design, c("2 candidates", "3 candidates", "4 candidates")
  )
  design_rank[is.na(design_rank)] <- 4L
  answer[order(answer$scenario, design_rank,
               answer$maximum_sample_size), ]
}

# A restrained, color-blind-friendly palette used consistently across figures.
combo_series_colors <- function() {
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
}

# Plot all six performance metrics against maximum sample size for one truth
# scenario. Design is encoded by both color and point shape. The separate
# legend band and large type keep the panels readable in the knitted HTML.
plot_sample_size_metrics <- function(summary_table, scenario_name) {
  x <- summary_table[summary_table$scenario == scenario_name, , drop = FALSE]
  if (!nrow(x)) stop("scenario_name is not present in summary_table")
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
  design_levels <- unique(x$design)
  colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")
  shapes <- c(16L, 17L, 15L, 18L)
  sizes <- sort(unique(x$maximum_sample_size))
  font_scale <- 1.4
  old <- par(no.readonly = TRUE)
  layout(
    matrix(c(1, 1, 2, 3, 4, 5, 6, 7), nrow = 4L, byrow = TRUE),
    heights = c(0.62, 1, 1, 1)
  )
  on.exit({
    layout(1)
    par(old)
  })

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center", legend = design_levels, col = colors,
    pch = shapes, lty = 1L, lwd = 3, pt.cex = 1.45,
    ncol = 4L, bty = "n", cex = font_scale * 1.16,
    title = scenario_name, title.adj = 0.5
  )

  for (k in seq_along(metric_names)) {
    values <- x[[metric_names[k]]]
    upper <- min(1, max(0.08, values, na.rm = TRUE) * 1.15)
    # Doubled labels need additional room so axis titles do not cover ticks.
    par(mar = c(7.2, 9.0, 4.4, 1.3), mgp = c(5.3, 1.6, 0))
    plot(
      NA, xlim = range(sizes), ylim = c(0, upper),
      xlab = "Maximum sample size", ylab = "Probability / proportion",
      main = metric_labels[k], xaxt = "n", yaxt = "n",
      cex.lab = font_scale * 1.36, cex.main = font_scale * 1.34
    )
    axis(1L, at = sizes, labels = sizes, cex.axis = font_scale * 1.18)
    axis(2L, las = 1L, cex.axis = font_scale * 1.18)
    abline(h = pretty(c(0, upper)), col = "#E5E7EB", lwd = 1)
    for (i in seq_along(design_levels)) {
      rows <- x$design == design_levels[i]
      one <- x[rows, , drop = FALSE]
      one <- one[order(one$maximum_sample_size), , drop = FALSE]
      lines(
        one$maximum_sample_size, one[[metric_names[k]]],
        type = "b", col = colors[i], pch = shapes[i],
        lwd = 2.8, cex = 1.35
      )
    }
  }
}

# Show the accuracy-safety trade-off as sample size increases. Each path runs
# from N=48 to N=90; moving upward improves PCS and moving left reduces
# allocation to combinations with true DLT probability above 0.33.
plot_sample_size_tradeoff <- function(summary_table) {
  scenario_levels <- unique(summary_table$scenario)
  design_levels <- unique(summary_table$design)
  colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")
  shapes <- c(16L, 17L, 15L, 18L)
  font_scale <- 1.4
  old <- par(no.readonly = TRUE)
  layout(
    matrix(c(1, 1, 2, 3, 4, 5), nrow = 3L, byrow = TRUE),
    heights = c(0.38, 1, 1)
  )
  on.exit({
    layout(1)
    par(old)
  })

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center", legend = design_levels, col = colors,
    pch = shapes, lty = 1L, lwd = 3, pt.cex = 1.45,
    ncol = 4L, bty = "n", cex = font_scale * 1.18,
    title = "Path from N = 48 to N = 90", title.adj = 0.5
  )

  for (scenario_name in scenario_levels) {
    x <- summary_table[
      summary_table$scenario == scenario_name, , drop = FALSE
    ]
    x_range <- range(x$B2_toxic_allocation, na.rm = TRUE)
    y_range <- range(x$PCS, na.rm = TRUE)
    x_pad <- max(0.012, diff(x_range) * 0.14)
    y_pad <- max(0.025, diff(y_range) * 0.14)
    x_limits <- c(max(0, x_range[1L] - x_pad),
                  min(1, x_range[2L] + x_pad))
    par(mar = c(5.5, 5.8, 3.8, 1.1))
    plot(
      NA, xlim = x_limits,
      ylim = c(max(0, y_range[1L] - y_pad),
               min(1, y_range[2L] + y_pad)),
      xlab = "B2 toxic allocation (lower is better)",
      ylab = "PCS (higher is better)", main = scenario_name,
      cex.lab = font_scale * 1.30, cex.main = font_scale * 1.30,
      cex.axis = font_scale * 1.16
    )
    grid(col = "#E5E7EB", lwd = 1)
    for (i in seq_along(design_levels)) {
      one <- x[x$design == design_levels[i], , drop = FALSE]
      one <- one[order(one$maximum_sample_size), , drop = FALSE]
      lines(
        one$B2_toxic_allocation, one$PCS,
        col = colors[i], lwd = 2.8
      )
      points(
        one$B2_toxic_allocation, one$PCS,
        pch = shapes[i], col = colors[i], cex = 1.35
      )
      text(
        one$B2_toxic_allocation[c(1L, nrow(one))],
        one$PCS[c(1L, nrow(one))],
        labels = c("48", "90"), col = colors[i],
        pos = c(2L, 4L), cex = font_scale * 0.90, font = 2L
      )
    }
  }
}

# Show a true two-drug toxicity surface on its dose grid. Percentages are
# printed directly in the cells; a smooth light line marks the target contour,
# and large circles identify the discrete combinations closest to the target.
plot_combo_scenario_grid <- function(probability, target, main) {
  n_a <- nrow(probability)
  n_b <- ncol(probability)
  plot(
    NA, xlim = c(0.5, n_b + 0.5), ylim = c(0.5, n_a + 0.5),
    xaxt = "n", yaxt = "n", xlab = "Drug B dose level",
    ylab = "Drug A dose level", main = main, asp = 1,
    cex.lab = 1.50, cex.main = 1.45
  )
  axis(1L, at = seq_len(n_b), labels = paste0("B", seq_len(n_b)),
       cex.axis = 1.35)
  axis(2L, at = seq_len(n_a), labels = paste0("A", seq_len(n_a)),
       las = 1L, cex.axis = 1.35)
  abline(v = seq(0.5, n_b + 0.5, by = 1),
         h = seq(0.5, n_a + 0.5, by = 1), col = "#D1D5DB", lwd = 1.2)

  # Interpolate the true 25% contour between the prespecified grid cells.
  contour(
    x = seq_len(n_b), y = seq_len(n_a), z = t(probability),
    levels = target, add = TRUE, drawlabels = FALSE,
    col = "#5BC0BE", lwd = 5
  )
  mtd <- true_combo_mtds(probability, target)
  points(
    mtd[, 2L], mtd[, 1L], pch = 21L, cex = 3.2,
    lwd = 2.2, bg = "#D5F5F2", col = "#168C8A"
  )
  grid <- expand.grid(drug_a = seq_len(n_a), drug_b = seq_len(n_b))
  values <- probability[cbind(grid$drug_a, grid$drug_b)]
  is_mtd <- paste(grid$drug_a, grid$drug_b) %in%
    paste(mtd[, 1L], mtd[, 2L])
  text(
    grid$drug_b, grid$drug_a, labels = paste0(round(100 * values), "%"),
    cex = 1.20, font = ifelse(is_mtd, 2L, 1L), col = "#111827"
  )
  legend("bottomleft", inset = 0.01, legend = "True 25% MTD contour",
         col = "#5BC0BE", lwd = 5, bty = "n", cex = 1.12)
  box(col = "#9CA3AF", lwd = 1.2)
}

# Bubble map for a 4 x 4 probability or allocation surface. Circle size and
# direct percentage labels carry the value; the pale fill is only supportive,
# so the plot remains readable without relying on a dark color scale.
plot_combo_bubble_surface <- function(probability, main, true_mtds = NULL,
                                      scale_max = 0.50, truth = NULL,
                                      target = NULL,
                                      toxic_threshold = NULL,
                                      show_safety_regions = FALSE) {
  n_a <- nrow(probability)
  n_b <- ncol(probability)
  plot(
    NA, xlim = c(0.5, n_b + 0.5), ylim = c(0.5, n_a + 0.5),
    xaxt = "n", yaxt = "n", xlab = "Drug B", ylab = "Drug A",
    main = main, asp = 1, cex.lab = 1.42, cex.main = 1.32
  )
  axis(1L, at = seq_len(n_b), labels = paste0("B", seq_len(n_b)),
       cex.axis = 1.32)
  axis(2L, at = seq_len(n_a), labels = paste0("A", seq_len(n_a)),
       cex.axis = 1.32, las = 1L)

  # In allocation plots, lightly mark the two safety regions used by the
  # operating characteristics: above-target cells for C1 and clearly toxic
  # cells for B2. C2 is an allocation-count metric, not a dose-grid region.
  if (show_safety_regions) {
    stopifnot(!is.null(truth), !is.null(target),
              !is.null(toxic_threshold))
    for (a in seq_len(n_a)) for (b in seq_len(n_b)) {
      fill <- if (truth[a, b] > toxic_threshold) {
        "#FDE2E2"
      } else if (truth[a, b] > target) {
        "#FFF1CC"
      } else {
        NA_character_
      }
      if (!is.na(fill)) {
        rect(b - 0.49, a - 0.49, b + 0.49, a + 0.49,
             col = fill, border = NA)
      }
    }
  }
  abline(v = seq(0.5, n_b + 0.5, by = 1),
         h = seq(0.5, n_a + 0.5, by = 1), col = "#E5E7EB", lwd = 1)
  grid <- expand.grid(drug_a = seq_len(n_a), drug_b = seq_len(n_b))
  values <- probability[cbind(grid$drug_a, grid$drug_b)]
  positive <- values >= 0.005
  point_size <- 1.45 + 4.3 * sqrt(pmin(values[positive] / scale_max, 1))
  points(
    grid$drug_b[positive], grid$drug_a[positive], pch = 21L,
    cex = point_size,
    bg = grDevices::adjustcolor("#56B4E9", alpha.f = 0.35),
    col = "#0072B2", lwd = 1.3
  )
  text(
    grid$drug_b[positive], grid$drug_a[positive],
    labels = paste0(round(100 * values[positive]), "%"),
    cex = 1.18, col = "#111827", font = 2L
  )
  if (!is.null(true_mtds)) {
    apply(true_mtds, 1L, function(cell) {
      rect(cell[2L] - 0.43, cell[1L] - 0.43,
           cell[2L] + 0.43, cell[1L] + 0.43,
           border = "#D55E00", lwd = 2.8)
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
  old <- par(no.readonly = TRUE)
  layout(
    matrix(c(1, 1, 2, 3, 4, 5, 6, 7), nrow = 4L, byrow = TRUE),
    heights = c(0.60, 1, 1, 1)
  )
  on.exit({
    layout(1)
    par(old)
  })

  # A dedicated legend row keeps the long scenario names out of the data area.
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center", legend = scenario_levels, col = colors,
    pch = shapes, lty = 1L, lwd = 3, pt.cex = 1.5,
    ncol = 2L, bty = "n", cex = 1.25
  )

  design_labels <- ifelse(
    grepl(" candidates$", design_levels),
    paste("Size", sub(" candidates$", "", design_levels)),
    sub("^Local parallel ", "Local\n", design_levels)
  )
  for (k in seq_along(metric_names)) {
    par(mar = c(5.8, 6.0, 3.7, 1.2))
    values <- summary_table[[metric_names[k]]]
    upper <- min(1, max(0.10, values, na.rm = TRUE) * 1.18)
    plot(
      NA, xlim = c(0.75, length(design_levels) + 0.25),
      ylim = c(0, upper), xaxt = "n", yaxt = "n", xlab = "Design",
      ylab = "Probability / proportion", main = metric_labels[k],
      cex.lab = 1.42, cex.main = 1.42
    )
    axis(1L, at = seq_along(design_levels),
         labels = design_labels, cex.axis = 1.28)
    axis(2L, las = 1L, cex.axis = 1.28)
    abline(h = pretty(c(0, upper)), col = "#E5E7EB", lwd = 1)
    for (i in seq_along(scenario_levels)) {
      row <- summary_table$scenario == scenario_levels[i]
      y_values <- summary_table[[metric_names[k]]][row]
      lines(seq_along(design_levels), y_values, type = "b",
            col = colors[i], pch = shapes[i], lwd = 3.0, cex = 1.80)
    }
  }
}

# Plot the four designs for one scenario in a two-column bubble-map page.
plot_combo_selection_bubbles <- function(results, scenarios, target,
                                         scenario_name) {
  design_levels <- unique(results$design)
  old <- par(
    mfrow = c(2, 2), mar = c(4.5, 4.8, 3.7, 1), oma = c(0, 0, 2.5, 0)
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
    outer = TRUE, cex = 1.48, font = 2L
  )
}

# Plot mean patient allocation for one scenario in a two-column page.
plot_combo_allocation_bubbles <- function(results, scenarios, scenario_name,
                                          target, toxic_threshold) {
  design_levels <- unique(results$design)
  allocation_names <- grep("^alloc_A[0-9]+_B[0-9]+$", names(results),
                           value = TRUE)
  old <- par(
    mfrow = c(2, 2), mar = c(4.5, 4.8, 3.7, 1), oma = c(0, 0, 2.5, 0)
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
      allocation_matrix, main = design_name, scale_max = 0.50,
      truth = truth, target = target, toxic_threshold = toxic_threshold,
      show_safety_regions = TRUE
    )
  }
  mtext(
    paste0(scenario_name,
           ": allocation; yellow = above target, pink = clearly toxic"),
    outer = TRUE, cex = 1.48, font = 2L
  )
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
  radius <- ifelse(visit_total > 1,
                   pmin(0.25, 0.14 + 0.012 * visit_total), 0)
  plot_a <- h$drug_a + radius * sin(angle)
  plot_b <- h$drug_b + radius * cos(angle)

  plot(
    plot_b, plot_a, type = "n",
    xlim = c(0.7, ncol(trial$true_probability) + 0.3),
    ylim = c(0.7, nrow(trial$true_probability) + 0.3),
    xlab = "Drug B level", ylab = "Drug A level", main = main,
    cex.lab = 1.50, cex.main = 1.45,
    xaxt = "n", yaxt = "n"
  )
  axis(1L, at = seq_len(ncol(trial$true_probability)),
       labels = paste0("B", seq_len(ncol(trial$true_probability))),
       cex.axis = 1.35)
  axis(2L, at = seq_len(nrow(trial$true_probability)),
       labels = paste0("A", seq_len(nrow(trial$true_probability))),
       cex.axis = 1.35, las = 1L)
  abline(v = seq_len(ncol(trial$true_probability)),
         h = seq_len(nrow(trial$true_probability)),
         col = "#EEF2F7", lwd = 1)

  # Display the true target contour behind the realized trial path.
  contour(
    x = seq_len(ncol(trial$true_probability)),
    y = seq_len(nrow(trial$true_probability)),
    z = t(trial$true_probability), levels = trial$target,
    add = TRUE, drawlabels = FALSE, col = "#A7E3CF", lwd = 7
  )
  true_mtds <- true_combo_mtds(trial$true_probability, trial$target)
  points(true_mtds[, 2L], true_mtds[, 1L], pch = 21L, cex = 3.5,
         bg = "#E2F7EF", col = "#65B891", lwd = 2.4)

  has_parallel_branches <- all(c(
    "parallel_block", "parallel_role",
    "parallel_parent_a", "parallel_parent_b"
  ) %in% names(h)) && any(grepl("^branch ", h$parallel_role))

  if (!has_parallel_branches) {
    lines(plot_b, plot_a, lwd = 1.8, col = "#64748B")
  } else {
    # The solid line follows the adaptive current path. Concurrent branch
    # cohorts are excluded so two sibling candidates are never drawn as a
    # sequential move from one to the other.
    path_rows <- which(!grepl("^branch ", h$parallel_role))
    if (length(path_rows) >= 2L) {
      lines(plot_b[path_rows], plot_a[path_rows],
            lwd = 1.8, col = "#64748B")
    }

    # Each parallel candidate receives its own dashed spoke from the common
    # launch combination. For example, escalation at A2B2 draws A2B2--A3B2
    # and A2B2--A2B3, never A3B2--A2B3.
    branch_rows <- which(grepl("^branch ", h$parallel_role))
    for (branch_row in branch_rows) {
      parent_row <- which(
        h$parallel_block == h$parallel_block[branch_row] &
          h$parallel_role == "launch"
      )[1L]
      if (!is.na(parent_row)) {
        segments(
          plot_b[parent_row], plot_a[parent_row],
          plot_b[branch_row], plot_a[branch_row],
          col = "#7C9CB8", lwd = 2.2, lty = 2L
        )
      }
    }
  }
  points(
    plot_b, plot_a, pch = 21L, cex = 2.0,
    bg = ifelse(h$cohort_dlt > 0, "#FDBA74", "#BAE6FD"),
    col = "#475569", lwd = 1.4
  )
  text(plot_b, plot_a, labels = h$cohort, cex = 1.00,
       col = "#111827", font = 2L)
  if (!anyNA(trial$selected)) {
    points(trial$selected[2L], trial$selected[1L],
           pch = 4L, cex = 2.8, lwd = 3, col = "#2563EB")
  }
}
