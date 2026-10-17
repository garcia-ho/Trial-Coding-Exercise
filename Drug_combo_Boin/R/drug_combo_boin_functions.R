# Functions used by Drug_combo_Boin.Rmd.
# Educational implementation only; not validated for clinical trial conduct.

# Calculate the optimal BOIN escalation and de-escalation boundaries.
# The usual defaults are phi1 = 0.6 * target and phi2 = 1.4 * target.
boin_boundaries <- function(target, phi1 = 0.6 * target,
                            phi2 = 1.4 * target) {
  stopifnot(
    target > 0, target < 1,
    phi1 > 0, phi1 < target,
    phi2 > target, phi2 < 1
  )
  lambda_e <- log((1 - phi1) / (1 - target)) /
    log(target * (1 - phi1) / (phi1 * (1 - target)))
  lambda_d <- log((1 - target) / (1 - phi2)) /
    log(phi2 * (1 - target) / (target * (1 - phi2)))
  c(escalate = lambda_e, deescalate = lambda_d)
}

# Beta-binomial posterior summaries for one combination. Jeffreys' Beta(0.5,
# 0.5) prior is used, matching the common BOIN combination specification.
combo_posterior_summary <- function(y, n, boundaries, target,
                                    prior_a = 0.5, prior_b = 0.5) {
  shape1 <- prior_a + y
  shape2 <- prior_b + n - y
  c(
    mean = shape1 / (shape1 + shape2),
    desirable = pbeta(boundaries["deescalate"], shape1, shape2) -
      pbeta(boundaries["escalate"], shape1, shape2),
    overdose = 1 - pbeta(target, shape1, shape2)
  )
}

# Return admissible row/column neighbors. Escalation increases exactly one drug
# by one level; de-escalation decreases exactly one drug by one level.
combo_neighbors <- function(current, n_a, n_b, direction, eliminated) {
  a <- current[1L]
  b <- current[2L]
  candidates <- if (direction == "escalate") {
    rbind(c(a + 1L, b), c(a, b + 1L))
  } else {
    rbind(c(a - 1L, b), c(a, b - 1L))
  }
  valid <- candidates[, 1L] >= 1L & candidates[, 1L] <= n_a &
    candidates[, 2L] >= 1L & candidates[, 2L] <= n_b
  candidates <- candidates[valid, , drop = FALSE]
  if (!nrow(candidates)) return(matrix(integer(), ncol = 2L))
  allowed <- !eliminated[candidates]
  candidates[allowed, , drop = FALSE]
}

# Choose the neighboring combination having the largest posterior probability
# of lying inside the BOIN stay interval. Random selection resolves exact ties.
choose_combo_neighbor <- function(candidates, y, n, boundaries, target) {
  if (!nrow(candidates)) return(NULL)
  scores <- apply(candidates, 1L, function(cell) {
    combo_posterior_summary(
      y[cell[1L], cell[2L]], n[cell[1L], cell[2L]], boundaries, target
    )["desirable"]
  })
  best <- which(abs(scores - max(scores)) < 1e-12)
  candidates[sample(best, 1L), ]
}

# Identify all combinations whose true DLT probability is closest to target.
# A combination trial may have more than one true MTD on an MTD contour.
true_combo_mtds <- function(true_probability, target) {
  distance <- abs(true_probability - target)
  which(abs(distance - min(distance)) < 1e-12, arr.ind = TRUE)
}

# Simulate one two-drug BOIN combination trial.
run_combo_boin_trial <- function(
    true_probability,
    target = 0.30,
    n_patients = 36L,
    cohort_size = 3L,
    start = c(1L, 1L),
    elimination_cutoff = 0.95,
    seed = NULL) {
  stopifnot(is.matrix(true_probability), all(true_probability > 0),
            all(true_probability < 1), length(start) == 2L)
  if (!is.null(seed)) set.seed(seed)

  n_a <- nrow(true_probability)
  n_b <- ncol(true_probability)
  boundaries <- boin_boundaries(target)
  n <- matrix(0L, n_a, n_b)
  y <- matrix(0L, n_a, n_b)
  eliminated <- matrix(FALSE, n_a, n_b)
  current <- as.integer(start)
  history <- data.frame()
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
    observed <- y[current[1L], current[2L]] /
      n[current[1L], current[2L]]
    posterior <- combo_posterior_summary(
      y[current[1L], current[2L]], n[current[1L], current[2L]],
      boundaries, target
    )

    # Eliminate the current and northeast combinations when excessive toxicity
    # is sufficiently probable. Eliminating (1,1) terminates the trial.
    if (n[current[1L], current[2L]] >= 3L &&
        posterior["overdose"] > elimination_cutoff) {
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
    if (eliminated[1L, 1L]) {
      action <- "stop: lowest combination too toxic"
      stopped_early <- TRUE
    } else if (action %in% c("escalate", "deescalate")) {
      candidates <- combo_neighbors(
        current, n_a, n_b, action, eliminated
      )
      chosen <- choose_combo_neighbor(candidates, y, n, boundaries, target)
      if (is.null(chosen)) {
        if (eliminated[current[1L], current[2L]]) {
          action <- "stop: no admissible lower neighbor"
          stopped_early <- TRUE
        } else {
          action <- "stay: boundary of dose grid"
        }
      } else {
        next_combo <- chosen
      }
    }

    history <- rbind(history, data.frame(
      cohort = cohort,
      first_patient = patient + 1L,
      last_patient = patient + cohort_n,
      drug_a = current[1L],
      drug_b = current[2L],
      cohort_dlt = sum(outcomes),
      cumulative_n = n[current[1L], current[2L]],
      cumulative_dlt = y[current[1L], current[2L]],
      observed_rate = observed,
      action = action,
      next_a = next_combo[1L],
      next_b = next_combo[2L]
    ))
    patient <- patient + cohort_n
    if (stopped_early) break
    current <- next_combo
  }

  posterior_mean <- (y + 0.5) / (n + 1)
  selectable <- n > 0 & !eliminated
  selected <- c(NA_integer_, NA_integer_)
  if (any(selectable)) {
    distance <- abs(posterior_mean - target)
    distance[!selectable] <- Inf
    best <- which(abs(distance - min(distance)) < 1e-12, arr.ind = TRUE)
    selected <- best[sample(seq_len(nrow(best)), 1L), ]
  }

  list(
    history = history,
    n = n,
    y = y,
    posterior_mean = posterior_mean,
    eliminated = eliminated,
    selected = unname(selected),
    boundaries = boundaries,
    target = target,
    true_probability = true_probability,
    true_mtds = true_combo_mtds(true_probability, target),
    stopped_early = stopped_early
  )
}

# Simulate the combination design under several true toxicity surfaces.
simulate_combo_boin_scenarios <- function(
    scenarios,
    target = 0.30,
    n_sim = 100L,
    n_patients = 36L,
    cohort_size = 3L,
    seed = 20261017L) {
  set.seed(seed)
  records <- vector("list", length(scenarios) * n_sim)
  index <- 0L

  for (scenario_name in names(scenarios)) {
    truth <- scenarios[[scenario_name]]
    true_mtds <- true_combo_mtds(truth, target)
    for (simulation in seq_len(n_sim)) {
      trial <- run_combo_boin_trial(
        truth, target, n_patients, cohort_size
      )
      h <- trial$history
      selected <- trial$selected
      correct <- !anyNA(selected) && any(
        true_mtds[, 1L] == selected[1L] &
          true_mtds[, 2L] == selected[2L]
      )
      true_at_treated <- truth[cbind(h$drug_a, h$drug_b)]
      mtd_treated <- vapply(seq_len(nrow(h)), function(i) {
        any(true_mtds[, 1L] == h$drug_a[i] &
              true_mtds[, 2L] == h$drug_b[i])
      }, logical(1L))

      index <- index + 1L
      records[[index]] <- data.frame(
        scenario = scenario_name,
        simulation = simulation,
        selected_a = selected[1L],
        selected_b = selected[2L],
        correct_selection = correct,
        dlt_rate = sum(h$cohort_dlt) / sum(h$last_patient - h$first_patient + 1L),
        mtd_allocation = sum(
          (h$last_patient - h$first_patient + 1L) * mtd_treated
        ) / sum(h$last_patient - h$first_patient + 1L),
        overdose_allocation = sum(
          (h$last_patient - h$first_patient + 1L) *
            (true_at_treated > target + 0.10)
        ) / sum(h$last_patient - h$first_patient + 1L),
        stopped_early = trial$stopped_early
      )
    }
  }
  do.call(rbind, records)
}

# Summarize operating characteristics by toxicity scenario.
summarize_combo_boin <- function(results) {
  groups <- split(results, results$scenario)
  answer <- do.call(rbind, lapply(groups, function(x) data.frame(
    scenario = x$scenario[1L],
    probability_correct_selection = mean(x$correct_selection),
    mean_dlt_rate = mean(x$dlt_rate),
    mean_mtd_allocation = mean(x$mtd_allocation),
    mean_overdose_allocation = mean(x$overdose_allocation),
    early_stop_probability = mean(x$stopped_early)
  )))
  rownames(answer) <- NULL
  answer
}

# Draw a labelled heatmap for a true or estimated toxicity surface.
plot_combo_surface <- function(probability, target, main,
                               selected = NULL, true_mtds = NULL) {
  n_a <- nrow(probability)
  n_b <- ncol(probability)
  colors <- hcl.colors(30L, "YlOrRd", rev = TRUE)
  image(
    seq_len(n_b), seq_len(n_a), t(probability),
    xlab = "Drug B dose level", ylab = "Drug A dose level",
    main = main, axes = FALSE, zlim = c(0, max(0.65, probability)),
    col = colors
  )
  axis(1L, at = seq_len(n_b))
  axis(2L, at = seq_len(n_a))
  labels <- formatC(as.vector(t(probability)), format = "f", digits = 2L)
  text(
    rep(seq_len(n_b), times = n_a),
    rep(seq_len(n_a), each = n_b),
    labels = labels, font = 2L
  )
  if (!is.null(true_mtds)) {
    apply(true_mtds, 1L, function(cell) {
      rect(cell[2L] - 0.45, cell[1L] - 0.45,
           cell[2L] + 0.45, cell[1L] + 0.45,
           border = "dodgerblue4", lwd = 3L)
    })
  }
  if (!is.null(selected) && !anyNA(selected)) {
    points(selected[2L], selected[1L], pch = 4L, cex = 2.2, lwd = 3L)
  }
  box()
}

# Plot the sequence of visited drug combinations for one simulated trial.
plot_combo_trial_path <- function(trial) {
  h <- trial$history
  plot(
    h$drug_b, h$drug_a, type = "o", lwd = 1.5, pch = 21L,
    bg = ifelse(h$cohort_dlt > 0, "firebrick", "white"),
    xlim = c(0.7, ncol(trial$true_probability) + 0.3),
    ylim = c(0.7, nrow(trial$true_probability) + 0.3),
    xlab = "Drug B dose level", ylab = "Drug A dose level",
    main = "Cohort allocation path", xaxt = "n", yaxt = "n"
  )
  axis(1L, at = seq_len(ncol(trial$true_probability)))
  axis(2L, at = seq_len(nrow(trial$true_probability)))
  text(h$drug_b, h$drug_a, labels = h$cohort, pos = 3L, cex = 0.8)
  legend(
    "topleft", c("No cohort DLT", "At least one cohort DLT"),
    pt.bg = c("white", "firebrick"), pch = 21L, bty = "n", cex = 0.8
  )
}

# Plot operating metrics for each simulation scenario.
plot_combo_operating_metrics <- function(summary_table) {
  metrics <- c(
    "probability_correct_selection",
    "mean_mtd_allocation",
    "mean_overdose_allocation",
    "mean_dlt_rate"
  )
  labels <- c("Correct selection", "MTD allocation",
              "Overdose allocation", "DLT rate")
  values <- t(as.matrix(summary_table[metrics]))
  colnames(values) <- summary_table$scenario
  barplot(
    values, beside = TRUE, ylim = c(0, 1),
    ylab = "Probability / proportion",
    main = "Combination-BOIN operating characteristics",
    legend.text = labels,
    args.legend = list(x = "topright", bty = "n", cex = 0.8)
  )
}

# Plot final combination-selection probabilities for all scenarios.
plot_combo_selection <- function(results, scenarios) {
  old <- par(mfrow = c(1, length(scenarios)), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  for (scenario_name in names(scenarios)) {
    truth <- scenarios[[scenario_name]]
    selected <- results[results$scenario == scenario_name &
                          !is.na(results$selected_a), ]
    probability <- matrix(0, nrow(truth), ncol(truth))
    if (nrow(selected)) {
      for (i in seq_len(nrow(selected))) {
        probability[selected$selected_a[i], selected$selected_b[i]] <-
          probability[selected$selected_a[i], selected$selected_b[i]] +
          1 / sum(results$scenario == scenario_name)
      }
    }
    plot_combo_surface(
      probability, target = 1, main = scenario_name,
      true_mtds = true_combo_mtds(truth, 0.30)
    )
  }
}
