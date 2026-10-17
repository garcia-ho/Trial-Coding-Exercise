# Functions used by CRM.Rmd.
# Educational implementation only; not validated for clinical trial conduct.

# Keep probabilities away from 0 and 1 before evaluating logarithms.
clamp_probability <- function(x, epsilon = 1e-12) {
  pmin(pmax(x, epsilon), 1 - epsilon)
}

# Construct a power-model skeleton from the Lee-Cheung indifference interval.
# target is phi, delta is the interval half-width, and prior_mtd is nu.
lee_cheung_skeleton <- function(target, delta, n_dose, prior_mtd) {
  stopifnot(
    target > 0, target < 1, delta > 0,
    target - delta > 0, target + delta < 1,
    n_dose >= 2L, prior_mtd >= 1L, prior_mtd <= n_dose
  )

  skeleton <- numeric(n_dose)
  skeleton[prior_mtd] <- target

  if (prior_mtd > 1L) {
    for (j in seq.int(prior_mtd, 2L)) {
      skeleton[j - 1L] <- exp(
        log(target - delta) * log(skeleton[j]) / log(target + delta)
      )
    }
  }
  if (prior_mtd < n_dose) {
    for (j in seq.int(prior_mtd, n_dose - 1L)) {
      skeleton[j + 1L] <- exp(
        log(target + delta) * log(skeleton[j]) / log(target - delta)
      )
    }
  }
  unname(skeleton)
}

# Evaluate the one-parameter power model p_j(alpha) = p_j ^ exp(alpha).
# Rows correspond to alpha values; columns correspond to dose levels.
crm_power_probability <- function(skeleton, alpha) {
  stopifnot(all(skeleton > 0), all(skeleton < 1))
  exp(outer(exp(alpha), log(skeleton), FUN = "*"))
}

# CRM log-likelihood, aggregated by dose for efficient repeated evaluation.
crm_log_likelihood <- function(alpha, dose, dlt, skeleton) {
  stopifnot(length(dose) == length(dlt), all(dlt %in% c(0, 1)))
  stopifnot(all(dose %in% seq_along(skeleton)))
  if (!length(dose)) return(rep(0, length(alpha)))

  n <- tabulate(dose, nbins = length(skeleton))
  y <- tabulate(dose[dlt == 1], nbins = length(skeleton))
  p <- clamp_probability(crm_power_probability(skeleton, alpha))
  rowSums(
    sweep(log(p), 2L, y, "*") +
      sweep(log1p(-p), 2L, n - y, "*")
  )
}

# Bayesian CRM fit by grid integration. Grid fitting is deterministic and fast,
# so it is used for interim decisions and operating-characteristic simulations.
crm_posterior_grid <- function(
    dose, dlt, skeleton, prior_mean = 0, prior_sd = sqrt(2),
    alpha_grid = seq(-4, 4, length.out = 801L)) {
  log_kernel <- crm_log_likelihood(alpha_grid, dose, dlt, skeleton) +
    dnorm(alpha_grid, prior_mean, prior_sd, log = TRUE)
  maximum <- max(log_kernel)
  unnormalized <- exp(log_kernel - maximum)
  mass <- unnormalized / sum(unnormalized)
  step <- alpha_grid[2L] - alpha_grid[1L]

  list(
    alpha = alpha_grid,
    density = mass / step,
    mass = mass,
    posterior_mean_alpha = sum(alpha_grid * mass),
    posterior_toxicity = colSums(crm_power_probability(skeleton, alpha_grid) * mass),
    log_marginal_likelihood = maximum + log(sum(unnormalized)) + log(step),
    skeleton = skeleton
  )
}

# Select the dose closest to target. During accrual, max_move = 1 prevents a
# jump of more than one level in either direction.
crm_recommend_dose <- function(posterior_toxicity, target, current_dose = NULL,
                               max_move = 1L) {
  desired <- which.min(abs(posterior_toxicity - target))
  if (is.null(current_dose) || is.infinite(max_move)) return(desired)
  lower <- max(1L, current_dose - max_move)
  upper <- min(length(posterior_toxicity), current_dose + max_move)
  min(max(desired, lower), upper)
}

# Generate independent binary DLT observations at one dose.
generate_dlt <- function(true_probability, n_patients) {
  stopifnot(true_probability >= 0, true_probability <= 1, n_patients >= 1L)
  rbinom(n_patients, 1L, true_probability)
}

# Build sequential cohort sizes. The first run_in_size patients are treated at
# the starting dose before the first model-based reassessment.
make_cohort_sizes <- function(n_patients, cohort_size, run_in_size) {
  stopifnot(
    n_patients >= 1L, cohort_size >= 1L, run_in_size >= 1L,
    run_in_size <= n_patients
  )
  sizes <- run_in_size
  remaining <- n_patients - run_in_size
  while (remaining > 0L) {
    next_size <- min(cohort_size, remaining)
    sizes <- c(sizes, next_size)
    remaining <- remaining - next_size
  }
  sizes
}

# Simulate one cohort-based, single-skeleton CRM trial.
run_crm_trial <- function(
    true_probabilities, skeleton, target, n_patients = 24L,
    cohort_size = 3L, start_dose = 1L, run_in_size = cohort_size,
    prior_sd = sqrt(2), seed = NULL) {
  stopifnot(length(true_probabilities) == length(skeleton))
  stopifnot(all(diff(true_probabilities) >= 0))
  if (!is.null(seed)) set.seed(seed)

  current <- start_dose
  history <- data.frame(
    patient = integer(), cohort = integer(), dose = integer(), dlt = integer()
  )
  cohort_sizes <- make_cohort_sizes(n_patients, cohort_size, run_in_size)
  starts <- cumsum(c(1L, head(cohort_sizes, -1L)))

  for (cohort in seq_along(cohort_sizes)) {
    first <- starts[cohort]
    cohort_n <- cohort_sizes[cohort]
    outcomes <- generate_dlt(true_probabilities[current], cohort_n)
    history <- rbind(history, data.frame(
      patient = seq.int(first, length.out = cohort_n), cohort = cohort,
      dose = current, dlt = outcomes
    ))
    fit <- crm_posterior_grid(history$dose, history$dlt, skeleton,
                              prior_sd = prior_sd)
    current <- crm_recommend_dose(fit$posterior_toxicity, target, current, 1L)
  }

  fit <- crm_posterior_grid(history$dose, history$dlt, skeleton,
                            prior_sd = prior_sd)
  list(
    history = history,
    fit = fit,
    selected_dose = crm_recommend_dose(fit$posterior_toxicity, target,
                                        max_move = Inf),
    allocation = tabulate(history$dose, nbins = length(skeleton)),
    true_mtd = which.min(abs(true_probabilities - target)),
    true_probabilities = true_probabilities,
    target = target
  )
}

# Random-walk Metropolis sampler for the CRM power parameter alpha.
# MCMC is used for learning and diagnostics; simulation decisions use the grid.
crm_mcmc <- function(
    dose, dlt, skeleton, n_iter = 12000L, burn = 2000L, thin = 5L,
    proposal_sd = 0.35, prior_mean = 0, prior_sd = sqrt(2), seed = NULL) {
  stopifnot(n_iter > burn, thin >= 1L, proposal_sd > 0)
  if (!is.null(seed)) set.seed(seed)

  log_posterior <- function(alpha) {
    crm_log_likelihood(alpha, dose, dlt, skeleton) +
      dnorm(alpha, prior_mean, prior_sd, log = TRUE)
  }

  chain <- numeric(n_iter)
  chain[1L] <- prior_mean
  accepted <- 0L
  current_log_posterior <- log_posterior(chain[1L])

  for (i in 2:n_iter) {
    proposal <- rnorm(1L, chain[i - 1L], proposal_sd)
    proposal_log_posterior <- log_posterior(proposal)
    if (log(runif(1L)) < proposal_log_posterior - current_log_posterior) {
      chain[i] <- proposal
      current_log_posterior <- proposal_log_posterior
      accepted <- accepted + 1L
    } else {
      chain[i] <- chain[i - 1L]
    }
  }

  samples <- chain[seq.int(burn + 1L, n_iter, by = thin)]
  list(
    samples = samples,
    full_chain = chain,
    acceptance_rate = accepted / (n_iter - 1L),
    posterior_toxicity = colMeans(crm_power_probability(skeleton, samples)),
    skeleton = skeleton
  )
}

# Fit all skeletons and average their posterior toxicity estimates using
# posterior model probabilities based on integrated marginal likelihoods.
bma_crm_posterior <- function(
    dose, dlt, skeletons, prior_model_weights = NULL,
    prior_sd = sqrt(2), alpha_grid = seq(-4, 4, length.out = 801L)) {
  stopifnot(length(skeletons) >= 2L)
  n_models <- length(skeletons)
  if (is.null(prior_model_weights)) {
    prior_model_weights <- rep(1 / n_models, n_models)
  }
  prior_model_weights <- prior_model_weights / sum(prior_model_weights)

  fits <- lapply(skeletons, function(skeleton) {
    crm_posterior_grid(dose, dlt, skeleton, prior_sd = prior_sd,
                       alpha_grid = alpha_grid)
  })
  log_weight <- log(prior_model_weights) +
    vapply(fits, function(x) x$log_marginal_likelihood, numeric(1L))
  log_weight <- log_weight - max(log_weight)
  model_weight <- exp(log_weight) / sum(exp(log_weight))
  toxicity <- vapply(
    fits, function(x) x$posterior_toxicity,
    numeric(length(skeletons[[1L]]))
  )

  list(
    fits = fits,
    posterior_model_weights = model_weight,
    posterior_toxicity_by_model = toxicity,
    posterior_toxicity = drop(toxicity %*% model_weight)
  )
}

# Simulate one BMA-CRM trial and retain model weights after every cohort.
run_bma_crm_trial <- function(
    true_probabilities, skeletons, target, n_patients = 24L,
    cohort_size = 3L, start_dose = 1L, run_in_size = cohort_size,
    prior_model_weights = NULL, prior_sd = sqrt(2), seed = NULL) {
  stopifnot(all(vapply(skeletons, length, integer(1L)) ==
                  length(true_probabilities)))
  if (!is.null(seed)) set.seed(seed)

  current <- start_dose
  history <- data.frame(
    patient = integer(), cohort = integer(), dose = integer(), dlt = integer()
  )
  cohort_sizes <- make_cohort_sizes(n_patients, cohort_size, run_in_size)
  starts <- cumsum(c(1L, head(cohort_sizes, -1L)))
  weight_history <- matrix(NA_real_, length(cohort_sizes), length(skeletons),
                           dimnames = list(NULL, names(skeletons)))

  for (cohort in seq_along(cohort_sizes)) {
    first <- starts[cohort]
    cohort_n <- cohort_sizes[cohort]
    outcomes <- generate_dlt(true_probabilities[current], cohort_n)
    history <- rbind(history, data.frame(
      patient = seq.int(first, length.out = cohort_n), cohort = cohort,
      dose = current, dlt = outcomes
    ))
    fit <- bma_crm_posterior(history$dose, history$dlt, skeletons,
                             prior_model_weights, prior_sd)
    weight_history[cohort, ] <- fit$posterior_model_weights
    current <- crm_recommend_dose(fit$posterior_toxicity, target, current, 1L)
  }

  fit <- bma_crm_posterior(history$dose, history$dlt, skeletons,
                           prior_model_weights, prior_sd)
  list(
    history = history,
    fit = fit,
    weight_history = weight_history,
    selected_dose = crm_recommend_dose(fit$posterior_toxicity, target,
                                        max_move = Inf),
    allocation = tabulate(history$dose, nbins = length(true_probabilities)),
    true_mtd = which.min(abs(true_probabilities - target)),
    true_probabilities = true_probabilities,
    target = target
  )
}

# Compare every single-skeleton CRM and BMA-CRM over multiple true scenarios.
simulate_crm_designs <- function(
    scenarios, skeletons, target, n_sim = 100L, n_patients = 24L,
    cohort_size = 3L, start_dose = 1L, seed = 20261004L,
    run_in_size = cohort_size) {
  set.seed(seed)
  design_names <- c(names(skeletons), "BMA-CRM")
  records <- vector("list", length(scenarios) * length(design_names) * n_sim)
  index <- 0L

  for (scenario_name in names(scenarios)) {
    truth <- scenarios[[scenario_name]]
    true_mtd <- which.min(abs(truth - target))
    for (design in design_names) {
      for (simulation in seq_len(n_sim)) {
        trial <- if (design == "BMA-CRM") {
          run_bma_crm_trial(truth, skeletons, target, n_patients,
                            cohort_size, start_dose,
                            run_in_size = run_in_size)
        } else {
          run_crm_trial(truth, skeletons[[design]], target, n_patients,
                        cohort_size, start_dose,
                        run_in_size = run_in_size)
        }
        index <- index + 1L
        records[[index]] <- data.frame(
          scenario = scenario_name, design = design, simulation = simulation,
          selected_dose = trial$selected_dose, true_mtd = true_mtd,
          dlt_rate = mean(trial$history$dlt),
          as.list(setNames(trial$allocation,
                           paste0("allocated_d", seq_along(trial$allocation)))),
          check.names = FALSE
        )
      }
    }
  }
  do.call(rbind, records)
}

# Cross three starting-dose run-in plans with all truth scenarios and
# skeletons. BMA is excluded so 3 plans x 3 scenarios x 3 skeletons produces
# the requested 27 design settings.
simulate_crm_run_in_plans <- function(
    scenarios, skeletons, target, run_in_sizes = c(3L, 6L, 9L),
    n_sim = 100L, n_patients = 24L, cohort_size = 3L,
    start_dose = 1L, seed = 20261004L) {
  stopifnot(length(run_in_sizes) == 3L)
  set.seed(seed)
  plan_names <- paste0(run_in_sizes, "-patient run-in")
  records <- vector(
    "list",
    length(run_in_sizes) * length(scenarios) * length(skeletons) * n_sim
  )
  index <- 0L

  for (plan_index in seq_along(run_in_sizes)) {
    run_in_size <- run_in_sizes[plan_index]
    for (scenario_name in names(scenarios)) {
      truth <- scenarios[[scenario_name]]
      true_mtd <- which.min(abs(truth - target))
      for (design in names(skeletons)) {
        for (simulation in seq_len(n_sim)) {
          trial <- run_crm_trial(
            truth, skeletons[[design]], target, n_patients,
            cohort_size, start_dose, run_in_size
          )
          index <- index + 1L
          records[[index]] <- data.frame(
            run_in_plan = plan_names[plan_index],
            run_in_size = run_in_size,
            scenario = scenario_name,
            design = design,
            simulation = simulation,
            selected_dose = trial$selected_dose,
            true_mtd = true_mtd,
            dlt_rate = mean(trial$history$dlt),
            as.list(setNames(
              trial$allocation,
              paste0("allocated_d", seq_along(trial$allocation))
            )),
            check.names = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, records)
}

# Summarize the 27 run-in x scenario x skeleton settings. Low-dose allocation
# is the proportion of patients assigned to dose levels 1 or 2.
summarize_crm_run_in_plans <- function(results, n_patients) {
  groups <- split(
    results,
    interaction(
      results$run_in_plan, results$scenario, results$design, drop = TRUE
    )
  )
  answer <- do.call(rbind, lapply(groups, function(x) data.frame(
    run_in_plan = x$run_in_plan[1L],
    run_in_size = x$run_in_size[1L],
    scenario = x$scenario[1L],
    design = x$design[1L],
    probability_correct_selection = mean(x$selected_dose == x$true_mtd),
    mean_dlt_rate = mean(x$dlt_rate),
    mean_low_dose_allocation = mean(
      (x$allocated_d1 + x$allocated_d2) / n_patients
    )
  )))
  rownames(answer) <- NULL
  answer
}

# Produce trial-level operating summaries for tables and plots.
summarize_crm_simulations <- function(results, n_dose, n_patients) {
  groups <- split(results, interaction(results$scenario, results$design,
                                       drop = TRUE))
  operating <- do.call(rbind, lapply(groups, function(x) data.frame(
    scenario = x$scenario[1L], design = x$design[1L],
    probability_correct_selection = mean(x$selected_dose == x$true_mtd),
    mean_dlt_rate = mean(x$dlt_rate)
  )))
  selection <- do.call(rbind, lapply(groups, function(x) data.frame(
    scenario = x$scenario[1L], design = x$design[1L], dose = seq_len(n_dose),
    probability = vapply(seq_len(n_dose),
                         function(d) mean(x$selected_dose == d), numeric(1L))
  )))
  allocation <- do.call(rbind, lapply(groups, function(x) data.frame(
    scenario = x$scenario[1L], design = x$design[1L], dose = seq_len(n_dose),
    proportion = colMeans(x[paste0("allocated_d", seq_len(n_dose))]) /
      n_patients
  )))
  rownames(operating) <- rownames(selection) <- rownames(allocation) <- NULL
  list(operating = operating, selection = selection, allocation = allocation)
}

# Visualize candidate skeletons.
plot_crm_skeletons <- function(skeletons, target) {
  values <- do.call(cbind, skeletons)
  matplot(seq_len(nrow(values)), values, type = "b", lty = 1L,
          pch = seq_len(ncol(values)), xlab = "Dose level",
          ylab = "Prior DLT probability", ylim = c(0, max(values, target)))
  abline(h = target, lty = 2L, col = "gray40")
  legend("topleft", names(skeletons), col = seq_along(skeletons),
         pch = seq_along(skeletons), lty = 1L, bty = "n")
}

# Visualize patient allocation; red points indicate DLTs.
plot_trial_path <- function(trial) {
  h <- trial$history
  plot(h$patient, h$dose, type = "s", lwd = 1.5,
       ylim = c(1, length(trial$true_probabilities)),
       xlab = "Patient", ylab = "Assigned dose", yaxt = "n")
  axis(2L, at = seq_along(trial$true_probabilities))
  points(h$patient, h$dose, pch = 19L,
         col = ifelse(h$dlt == 1L, "firebrick", "gray25"))
  legend("topleft", c("No DLT", "DLT"),
         col = c("gray25", "firebrick"), pch = 19L, bty = "n")
}

# Compare estimated and true DLT probabilities, with an optional additional
# curve such as a prior mean or a posterior estimate from another CRM model.
plot_posterior_toxicity <- function(
    estimate, truth, target, main = "Posterior toxicity",
    estimate_label = "Posterior mean", comparison = NULL,
    comparison_label = "Comparison") {
  dose <- seq_along(truth)
  plot(dose, truth, type = "b", pch = 19L,
       ylim = c(0, max(truth, estimate, comparison, target)), xlab = "Dose level",
       ylab = "DLT probability", main = main)
  lines(dose, estimate, type = "b", pch = 1L, col = "navy")
  if (!is.null(comparison)) {
    stopifnot(length(comparison) == length(truth))
    lines(dose, comparison, type = "b", pch = 2L, lty = 3L,
          col = "darkorange3")
  }
  abline(h = target, lty = 2L, col = "gray40")
  labels <- c("Truth", estimate_label)
  colors <- c("black", "navy")
  points <- c(19L, 1L)
  line_types <- c(1L, 1L)
  if (!is.null(comparison)) {
    labels <- c(labels, comparison_label)
    colors <- c(colors, "darkorange3")
    points <- c(points, 2L)
    line_types <- c(line_types, 3L)
  }
  legend(
    "topleft", c(labels, "Target"),
    col = c(colors, "gray40"), pch = c(points, NA),
    lty = c(line_types, 2L), bty = "n", cex = 0.85
  )
}

# Basic MCMC trace and marginal-posterior diagnostics.
plot_crm_mcmc <- function(fit) {
  old <- par(mfrow = c(1, 2)); on.exit(par(old))
  plot(fit$full_chain, type = "l", xlab = "Iteration",
       ylab = expression(alpha), main = "Metropolis trace")
  hist(fit$samples, breaks = 30L, probability = TRUE, col = "gray85",
       border = "white", xlab = expression(alpha), main = "Posterior samples")
  lines(density(fit$samples), col = "navy", lwd = 2L)
}

# Show how BMA posterior model weights change after each cohort.
plot_bma_weights <- function(trial) {
  matplot(seq_len(nrow(trial$weight_history)), trial$weight_history,
          type = "b", lty = 1L, pch = seq_len(ncol(trial$weight_history)),
          ylim = c(0, 1), xlab = "Cohort",
          ylab = "Posterior model probability")
  legend("topright", colnames(trial$weight_history),
         col = seq_len(ncol(trial$weight_history)),
         pch = seq_len(ncol(trial$weight_history)), lty = 1L, bty = "n")
}

# Plot correct-selection probability and average observed DLT rate.
plot_operating_metrics <- function(summary_object) {
  x <- summary_object$operating
  scenarios <- unique(x$scenario); designs <- unique(x$design)
  make_matrix <- function(variable) {
    out <- matrix(NA_real_, length(designs), length(scenarios),
                  dimnames = list(designs, scenarios))
    for (i in seq_len(nrow(x))) out[x$design[i], x$scenario[i]] <- x[[variable]][i]
    out
  }
  old <- par(mfrow = c(1, 2), mar = c(7, 4, 3, 1)); on.exit(par(old))
  barplot(make_matrix("probability_correct_selection"), beside = TRUE,
          ylim = c(0, 1), las = 2L, ylab = "Probability",
          main = "Correct MTD selection", legend.text = designs,
          args.legend = list(x = "topright", bty = "n", cex = 0.7))
  barplot(make_matrix("mean_dlt_rate"), beside = TRUE,
          ylim = c(0, max(0.5, x$mean_dlt_rate)), las = 2L,
          ylab = "Mean DLT rate", main = "Observed toxicity")
}

# Internal helper for grouped selection/allocation bar plots.
plot_dose_metric <- function(data, scenario, value, ylab, main) {
  x <- data[data$scenario == scenario, ]
  designs <- unique(x$design); doses <- sort(unique(x$dose))
  values <- matrix(0, length(designs), length(doses),
                   dimnames = list(designs, paste0("Dose ", doses)))
  for (i in seq_len(nrow(x))) {
    values[x$design[i], paste0("Dose ", x$dose[i])] <- x[[value]][i]
  }
  barplot(values, beside = TRUE, ylim = c(0, 1), xlab = "Dose",
          ylab = ylab, main = main, legend.text = designs,
          args.legend = list(x = "topright", bty = "n", cex = 0.7))
}

plot_selection_probabilities <- function(summary_object, scenario) {
  plot_dose_metric(summary_object$selection, scenario, "probability",
                   "Selection probability", scenario)
}

plot_allocation_proportions <- function(summary_object, scenario) {
  plot_dose_metric(summary_object$allocation, scenario, "proportion",
                   "Mean allocation proportion", scenario)
}

# Display the 3 x 3 x 3 factorial results compactly. Each panel is one
# skeleton; rows are truth scenarios and columns are starting-dose run-in plans.
plot_crm_factorial_heatmaps <- function(
    factorial_summary, metric, title, zlim = c(0, 1), digits = 2L) {
  designs <- unique(factorial_summary$design)
  scenarios <- unique(factorial_summary$scenario)
  plans <- unique(
    factorial_summary[order(factorial_summary$run_in_size), "run_in_plan"]
  )
  colors <- hcl.colors(20L, palette = "YlOrRd", rev = TRUE)
  old <- par(
    mfrow = c(1, length(designs)), mar = c(7, 7.5, 4.5, 1),
    oma = c(0, 0, 2.5, 0)
  )
  on.exit(par(old))

  for (design in designs) {
    values <- matrix(
      NA_real_, nrow = length(plans), ncol = length(scenarios),
      dimnames = list(plans, scenarios)
    )
    x <- factorial_summary[factorial_summary$design == design, ]
    for (i in seq_len(nrow(x))) {
      values[x$run_in_plan[i], x$scenario[i]] <- x[[metric]][i]
    }

    image(
      seq_along(plans), seq_along(scenarios), values,
      axes = FALSE, xlab = "", ylab = "",
      main = design, zlim = zlim, col = colors, cex.main = 1.3
    )
    axis(1L, at = seq_along(plans), labels = plans, las = 2L, cex.axis = 1.05)
    axis(2L, at = seq_along(scenarios), labels = scenarios,
         las = 2L, cex.axis = 1.05)
    text(
      rep(seq_along(plans), times = length(scenarios)),
      rep(seq_along(scenarios), each = length(plans)),
      labels = formatC(as.vector(values), format = "f", digits = digits),
      font = 2L, cex = 1.2
    )
    box()
  }
  mtext(title, outer = TRUE, line = 0.7, font = 2L, cex = 1.3)
}
