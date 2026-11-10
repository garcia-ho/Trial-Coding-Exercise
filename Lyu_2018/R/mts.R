# Appendix p. 18: sequence-specific TPC, possible orders, PAVA, then voting.
# Table 2 does not report its final MTS. This module demonstrates the
# published algorithm with its underspecified choices made explicit.

pava <- function(y, weights = rep(1, length(y))) {
  stopifnot(length(y) > 0, length(y) == length(weights),
            all(is.finite(y)), all(is.finite(weights)), all(weights > 0))
  value <- y
  weight <- weights
  first <- last <- seq_along(y)
  repeat {
    bad <- which(diff(value) < 0)
    if (!length(bad)) break
    b <- bad[1]
    # Merge the first adjacent violation, then recheck all blocks.
    value[b] <- sum(value[b:(b + 1)] * weight[b:(b + 1)]) /
                sum(weight[b:(b + 1)])
    weight[b] <- weight[b] + weight[b + 1]
    last[b] <- last[b + 1]
    value <- value[-(b + 1)]
    weight <- weight[-(b + 1)]
    first <- first[-(b + 1)]
    last <- last[-(b + 1)]
  }
  fitted <- numeric(length(y))
  for (b in seq_along(value)) fitted[first[b]:last[b]] <- value[b]
  fitted
}

permutations <- function(v) {
  if (length(v) == 1) return(matrix(v, nrow = 1))
  do.call(rbind, lapply(seq_along(v), function(i) {
    cbind(v[i], permutations(v[-i]))
  }))
}

sequence_summary <- function(fit, screen_min_n = fit$settings$safety_min_n) {
  stopifnot(length(screen_min_n) == 1, is.finite(screen_min_n),
            screen_min_n >= 0, screen_min_n %% 1 == 0)
  o <- fit$observed
  complete <- which(tabulate(o$patient, nbins = fit$n_patient) == fit$n_cycle)
  if (!length(complete)) return(data.frame())
  # Only patients completing all cycles contribute to sequence-specific TPC.
  doses <- matrix(NA_integer_, length(complete), fit$n_cycle)
  tox <- matrix(NA_integer_, length(complete), fit$n_cycle)
  for (i in seq_along(complete)) {
    r <- o[o$patient == complete[i], ]
    doses[i, r$cycle] <- r$dose
    tox[i, r$cycle] <- r$tox
  }
  key <- apply(doses, 1, paste, collapse = ",")
  rows <- lapply(unique(key), function(s) {
    idx <- which(key == s)
    path <- doses[idx[1], ]
    # The final posterior screen can have a different minimum n from the
    # intratrial Figure 2 override. Simulations pass this choice explicitly.
    unsafe <- any(screen_min_n <= fit$counts$n[cbind(path, seq_along(path))] &
      overdose_probability(fit$counts$x[cbind(path, seq_along(path))],
                           fit$counts$n[cbind(path, seq_along(path))],
                           fit$settings$target) > fit$settings$xi)
    data.frame(sequence = s, intensity = sum(path), patients = length(idx),
               events = sum(tox[idx, , drop = FALSE]),
               patient_cycles = length(idx) * fit$n_cycle,
               raw_tpc = mean(tox[idx, , drop = FALSE]),
               nonincreasing = all(diff(path) <= 0),
               admissible = !unsafe && !any(fit$snapshots[[length(fit$snapshots)]]$excluded[
                 cbind(path, seq_along(path))]))
  })
  tab <- do.call(rbind, rows)
  tab[order(tab$intensity, tab$sequence), ]
}

# weights and nearest_tie are REQUIRED: the appendix does not fully specify
# PAVA weights or nearest-target ties WITHIN a possible simple ordering.
select_mts <- function(fit, weights, nearest_tie, max_orders = 100000L,
                        screen_min_n = fit$settings$safety_min_n,
                        no_mts_fallback = "restricted_revote") {
  weights <- match.arg(weights, c("equal", "patient_cycles"))
  nearest_tie <- match.arg(nearest_tie, c("lower", "higher"))
  no_mts_fallback <- match.arg(no_mts_fallback, c("restricted_revote", "none"))
  tab <- sequence_summary(fit, screen_min_n = screen_min_n)
  if (nrow(tab)) tab <- tab[tab$admissible, ]
  empty <- list(recommendation = NA_character_, sequences = tab,
                votes = data.frame(), per_order = data.frame(),
                fitted = matrix(numeric(), 0, 0), fallback_used = FALSE)
  if (fit$stopped || !nrow(tab)) return(empty)
  rownames(tab) <- NULL
  groups <- split(seq_len(nrow(tab)), tab$intensity)
  n_orders <- prod(vapply(groups, function(g) factorial(length(g)), numeric(1)))
  if (n_orders > max_orders) stop("Too many possible orders; exact enumeration exceeds max_orders.")
  orders <- matrix(integer(), nrow = 1)
  for (g in groups) {
    perm <- permutations(g)
    orders <- do.call(rbind, lapply(seq_len(nrow(orders)), function(r) {
      cbind(matrix(rep(orders[r, ], each = nrow(perm)), nrow = nrow(perm)), perm)
    }))
  }
  w <- if (weights == "equal") rep(1, nrow(tab)) else tab$patient_cycles
  fitted <- matrix(NA_real_, nrow(orders), nrow(tab), dimnames = list(NULL, tab$sequence))
  chosen <- integer(nrow(orders))
  # Resolve nearest-target ties by position in this ordering. The simulation
  # uses the original higher-position rule, also during fallback revoting.
  nearest <- function(q, candidates = seq_along(q)) {
    distance <- abs(q[candidates] - fit$settings$target)
    tied <- candidates[abs(distance - min(distance)) < 1e-12]
    if (nearest_tie == "lower") min(tied) else max(tied)
  }
  for (r in seq_len(nrow(orders))) {
    ord <- orders[r, ]
    q <- pava(tab$raw_tpc[ord], w[ord])
    fitted[r, ord] <- q
    chosen[r] <- ord[nearest(q)]
  }
  votes <- tab[, c("sequence", "intensity", "nonincreasing")]
  votes$frequency <- tabulate(chosen, nbins = nrow(tab))
  # Appendix filters to nonincreasing paths AFTER obtaining per-order votes.
  eligible <- which(votes$nonincreasing & votes$frequency > 0)
  frequency <- votes$frequency
  votes$fallback_frequency <- 0L
  fallback_chosen <- rep(NA_integer_, nrow(orders))
  fallback_used <- FALSE
  # NO-MTS fallback (an explicit completion, not specified in the appendix):
  # Only if every original winner is increasing, repeat the nearest-target
  # vote within the observed, safe, nonincreasing set. Reuse the ORIGINAL
  # PAVA fits/orderings; do not refit, relax safety, or invent an unused path.
  if (!length(eligible) && any(votes$nonincreasing) &&
      no_mts_fallback == "restricted_revote") {
    fallback_used <- TRUE
    for (r in seq_len(nrow(orders))) {
      ord <- orders[r, ]
      positions <- which(votes$nonincreasing[ord])
      fallback_chosen[r] <- ord[nearest(fitted[r, ord], positions)]
    }
    votes$fallback_frequency <- tabulate(fallback_chosen, nbins = nrow(tab))
    frequency <- votes$fallback_frequency
    eligible <- which(votes$nonincreasing & frequency > 0)
  }
  recommendation <- NA_character_
  if (length(eligible)) {
    best <- eligible[frequency[eligible] == max(frequency[eligible])]
    # Final frequency ties: lowest intensity, then lowest first, second, ...
    # dose. Numeric lexicographic ordering avoids multi-digit string ordering.
    path <- do.call(rbind, strsplit(votes$sequence[best], ",", fixed = TRUE))
    storage.mode(path) <- "integer"
    tie_order <- do.call(order, c(list(votes$intensity[best]), as.data.frame(path)))
    recommendation <- votes$sequence[best[tie_order[1]]]
  }
  list(recommendation = recommendation, sequences = tab, votes = votes,
       per_order = data.frame(order = seq_len(nrow(orders)),
                              ordering = apply(orders, 1, function(o) paste(tab$sequence[o], collapse = " < ")),
                              selected = tab$sequence[chosen],
                              fallback_selected = tab$sequence[fallback_chosen]),
       fitted = fitted, weights = weights, nearest_tie = nearest_tie,
       no_mts_fallback = no_mts_fallback, fallback_used = fallback_used)
}
