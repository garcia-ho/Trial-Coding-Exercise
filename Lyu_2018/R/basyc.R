# BaSyc allocation, Appendix pp. 16-17.
# Keep observed data, cumulative dose/cycle counts, and decisions separate.

dose_cycle_counts <- function(observed, n_dose, n_cycle) {
  x <- n <- matrix(0L, n_dose, n_cycle,
                   dimnames = list(paste0("D", seq_len(n_dose)),
                                   paste0("C", seq_len(n_cycle))))
  for (r in seq_len(nrow(observed))) {
    d <- observed$dose[r]
    j <- observed$cycle[r]
    n[d, j] <- n[d, j] + 1L
    x[d, j] <- x[d, j] + observed$tox[r]
  }
  list(x = x, n = n)
}

# Exclusion is permanent and extends to higher doses AND later cycles.
update_safety <- function(counts, excluded, settings) {
  risk <- overdose_probability(counts$x, counts$n, settings$target)
  dim(risk) <- dim(counts$n)
  bad <- risk > settings$xi & counts$n >= settings$safety_min_n
  for (r in seq_len(nrow(bad))) {
    for (j in seq_len(ncol(bad))) {
      if (bad[r, j]) excluded[r:nrow(bad), j:ncol(bad)] <- TRUE
    }
  }
  list(excluded = excluded, stop_trial = any(bad[1, ]), risk = risk)
}

# DCR uses the LOWEST dose in the NEXT cycle, not the current patient's dose.
dcr <- function(cycle, counts, settings) {
  if (cycle == ncol(counts$n)) return("Exit")
  n <- counts$n[1, cycle + 1]
  risk <- overdose_probability(counts$x[1, cycle + 1], n, settings$target)
  if (n >= settings$safety_min_n && risk > settings$xi) "NG" else "G"
}

# One reusable mTPI-2 module, evaluated at the requested dose/cycle cell.
# decision_fun may be replaced by another verified mTPI-2 implementation
# with signature function(x, n, settings), returning E/S/D/DU.
dose_candidate <- function(dose, cycle, counts, settings,
                           decision_fun = mtpi2_decision) {
  x <- counts$x[dose, cycle]
  n <- counts$n[dose, cycle]
  action <- decision_fun(x, n, settings)
  stopifnot(action %in% c("E", "S", "D", "DU"))
  move <- switch(action, E = 1L, S = 0L, D = -1L, DU = -1L)
  candidate <- max(1L, min(settings$n_dose, dose + move))
  list(dose = candidate, action = action, x = x, n = n)
}

# Apply exclusion after the DAR minimum; never replace a forbidden dose
# by a higher dose. NA means no available dose at or below the candidate.
allowed_dose <- function(candidate, cycle, excluded) {
  available <- which(!excluded[, cycle] & seq_len(nrow(excluded)) <= candidate)
  if (length(available)) max(available) else NA_integer_
}

dar1 <- function(previous_dose, counts, excluded, settings,
                 decision_fun = mtpi2_decision) {
  z <- dose_candidate(previous_dose, 1, counts, settings, decision_fun)
  list(dose = allowed_dose(z$dose, 1, excluded), action = z$action,
       i1 = z$dose, i2 = NA_integer_,
       current = sprintf("D%d C1: %d/%d %s", previous_dose, z$x, z$n, z$action),
       previous = "")
}

dar2 <- function(current_dose, current_tox, cycle, previous_doses,
                 counts, excluded, settings, decision_fun = mtpi2_decision) {
  # If the previous cohort split, use ALL its dose levels in cycle j+1,
  # not just its most recent patient or one of its subcohorts.
  stopifnot(length(previous_doses) > 0)
  a <- dose_candidate(current_dose, cycle, counts, settings, decision_fun)
  b <- lapply(sort(unique(previous_doses)), dose_candidate, cycle = cycle + 1,
              counts = counts, settings = settings, decision_fun = decision_fun)
  i2 <- min(vapply(b, function(z) z$dose, integer(1)))
  # DAR-2a: min(i1,i2); DAR-2b: min(max(current dose-1,1),i2).
  personal <- if (current_tox == 0) a$dose else max(current_dose - 1L, 1L)
  result <- allowed_dose(min(personal, i2), cycle + 1, excluded)
  details <- vapply(seq_along(b), function(r) {
    sprintf("D%d C%d: %d/%d %s -> D%d", sort(unique(previous_doses))[r],
            cycle + 1, b[[r]]$x, b[[r]]$n, b[[r]]$action, b[[r]]$dose)
  }, character(1))
  list(dose = result, action = if (current_tox == 0) "DAR-2a" else "DAR-2b",
       i1 = a$dose, i2 = i2,
       current = sprintf("D%d C%d: %d/%d %s", current_dose, cycle, a$x, a$n, a$action),
       previous = paste(details, collapse = "; "))
}

# Replay binary outcomes, or simulate outcomes from pre-generated uniforms.
# With toxicity_prob supplied, outcomes[i,j] is U(0,1) and the observed DLT
# is I(U < p[d,j]) at the dose ACTUALLY assigned. Future outcomes stay hidden.
run_basyc <- function(outcomes, cohort_size = 3L, cycle_weeks = 8L,
                      settings = basyc_settings(), decision_fun = mtpi2_decision,
                      toxicity_prob = NULL, global_monitor = "continuous") {
  global_monitor <- match.arg(global_monitor, c("continuous", "after_first_sequence"))
  outcomes <- as.matrix(outcomes)
  stopifnot(nrow(outcomes) > 0, ncol(outcomes) > 0,
            cohort_size >= 1, cohort_size %% 1 == 0,
            nrow(outcomes) %% cohort_size == 0, cycle_weeks > 0)
  n_patient <- nrow(outcomes)
  n_cycle <- ncol(outcomes)
  n_cohort <- n_patient / cohort_size
  if (is.null(toxicity_prob)) {
    stopifnot(all(outcomes %in% c(0, 1)))
  } else {
    stopifnot(is.matrix(toxicity_prob),
              identical(dim(toxicity_prob), c(as.integer(settings$n_dose), n_cycle)),
              all(is.finite(toxicity_prob)), all(toxicity_prob >= 0 & toxicity_prob <= 1),
              all(is.finite(outcomes)), all(outcomes >= 0 & outcomes < 1))
  }
  cohort <- rep(seq_len(n_cohort), each = cohort_size)
  assignments <- data.frame(week = numeric(), patient = integer(), cohort = integer(),
                            cycle = integer(), dose = integer(), end_week = numeric(),
                            rule = character(), i1 = integer(), i2 = integer(),
                            current = character(), previous = character())
  observed <- data.frame(patient = integer(), cohort = integer(), cycle = integer(),
                         dose = integer(), end_week = numeric(), tox = integer())
  excluded <- matrix(FALSE, settings$n_dose, n_cycle)
  snapshots <- list()
  continuation <- data.frame(week = numeric(), patient = integer(),
                              cycle = integer(), decision = character())
  stopped <- FALSE
  global_start_week <- NA_real_
  stop_event <- list(cycle=NA_integer_,x=NA_integer_,n=NA_integer_,risk=NA_real_)
  counts <- dose_cycle_counts(observed, settings$n_dose, n_cycle)

  # Each epoch: observe ALL outcomes due now, update counts/safety, then
  # make decisions. No future-cycle outcome is available to a decision.
  for (epoch in 0:(n_cohort + n_cycle - 1L)) {
    week <- epoch * cycle_weeks
    due <- assignments[assignments$end_week == week, ]
    if (nrow(due)) {
      new <- due[, c("patient", "cohort", "cycle", "dose", "end_week")]
      new$tox <- outcomes[cbind(new$patient, new$cycle)]
      if (!is.null(toxicity_prob)) {
        new$tox <- as.integer(new$tox < toxicity_prob[cbind(new$dose, new$cycle)])
      }
      observed <- rbind(observed, new)
      counts <- dose_cycle_counts(observed, settings$n_dose, n_cycle)
    }

    # Experimental mode starts BOTH global termination and permanent exclusion
    # after cohort 1 completes every cycle. DCR remains active independently.
    first_sequence_complete <- sum(observed$cohort==1 & observed$cycle==n_cycle)==cohort_size
    global_ready <- global_monitor=="continuous" || first_sequence_complete
    more_treatment <- any(due$cycle<n_cycle) || epoch+1L<=n_cohort
    if (nrow(due) && global_ready) {
      if (is.na(global_start_week)) global_start_week <- week
      safety <- update_safety(counts,excluded,settings)
      excluded <- safety$excluded
      # In the experimental mode, a terminal visit cannot be an EARLY stop.
      # Its safety data still inform exclusion and the final MTS screen.
      stopped <- safety$stop_trial && (global_monitor=="continuous" || more_treatment)
      if (stopped) {
        j_stop <- unname(which(counts$n[1,]>=settings$safety_min_n &
                          safety$risk[1,]>settings$xi)[1])
        stop_event <- list(cycle=j_stop,x=counts$x[1,j_stop],n=counts$n[1,j_stop],
                           risk=safety$risk[1,j_stop])
      }
    }
    allocation_decision <- if (global_ready) decision_fun else function(x,n,s) {
      action <- decision_fun(x,n,s)
      if (action=="DU") "D" else action  # Exclusion has not started yet.
    }

    # Continue previously enrolled patients, subject to DCR and safety.
    if (nrow(due)) {
      for (r in seq_len(nrow(due))) {
        row <- due[r, ]
        j <- row$cycle
        gate <- if (stopped) "Stop trial" else dcr(j, counts, settings)
        continuation <- rbind(continuation,
                              data.frame(week, patient = row$patient,
                                         cycle = j, decision = gate))
        if (gate != "G") next
        if (row$cohort == 1) {
          # Appendix default: cohort 1 remains at D1 in every cycle.
          z <- list(dose = allowed_dose(1, j + 1, excluded), action = "Default",
                    i1 = NA_integer_, i2 = NA_integer_, current = "", previous = "")
        } else {
          prev <- observed[observed$cohort == row$cohort - 1 &
                             observed$cycle == j + 1, ]
          # Under equal eight-week cycles, the previous cohort finishes
          # the destination cycle at this epoch before decisions are made.
          if (!nrow(prev)) stop("Previous cohort has not completed the next cycle.")
          current_tox <- new$tox[new$patient == row$patient & new$cycle == j]
          z <- dar2(row$dose, current_tox, j, prev$dose,
                    counts, excluded, settings, allocation_decision)
        }
        if (is.na(z$dose)) stop("No available next-cycle dose after safety checks.")
        assignments <- rbind(assignments,
          data.frame(week, patient = row$patient, cohort = row$cohort,
                     cycle = j + 1L, dose = z$dose, end_week = week + cycle_weeks,
                     rule = z$action, i1 = z$i1, i2 = z$i2,
                     current = z$current, previous = z$previous))
      }
    }

    # Fast enrollment: a new cohort starts when its predecessor finishes C1.
    # Cohort 1 starts at week 0; the last cohort starts at week 32 in Table 2.
    k <- epoch + 1L
    if (!stopped && k <= n_cohort) {
      if (k == 1) {
        z <- list(dose = 1L, action = "Default", i1 = NA_integer_,
                  i2 = NA_integer_, current = "", previous = "")
      } else {
        last <- due[due$cohort == k - 1 & due$cycle == 1, ]
        stopifnot(length(unique(last$dose)) == 1)
        z <- dar1(last$dose[1], counts, excluded, settings, allocation_decision)
      }
      if (is.na(z$dose)) stop("No available cycle-1 dose after safety checks.")
      ids <- which(cohort == k)
      assignments <- rbind(assignments,
        data.frame(week, patient = ids, cohort = k, cycle = 1L, dose = z$dose,
                   end_week = week + cycle_weeks,
                   rule = if (k == 1) "Default" else "DAR-1",
                   i1 = z$i1, i2 = z$i2, current = z$current, previous = z$previous))
    }
    snapshots[[as.character(week)]] <- list(x = counts$x, n = counts$n,
                                            excluded = excluded,global_active=global_ready)
    if (stopped) break
  }
  rownames(assignments) <- rownames(observed) <- rownames(continuation) <- NULL
  list(assignments = assignments, observed = observed, counts = counts,
       continuation = continuation, snapshots = snapshots, stopped = stopped,
       end_week = week, settings = settings, n_cycle = n_cycle,
       n_patient = n_patient, cohort_size = cohort_size,
       global_monitor=global_monitor,global_start_week=global_start_week,
       stop_event=stop_event)
}
