# The small mTPI-2 module used by every BaSyc dose/cycle decision.
# Only base R is needed. This is mTPI-2, not an empirical x/n rule or BOIN.

basyc_settings <- function(n_dose = 3L, target = 0.30,
                           ei = c(0.25, 0.35), xi = 0.95,
                           safety_min_n = 0L) {
  stopifnot(length(n_dose) == 1, n_dose >= 1, n_dose %% 1 == 0,
            length(ei) == 2, 0 < ei[1], ei[1] < target,
            target < ei[2], ei[2] < 1, 0 < xi, xi < 1,
            safety_min_n >= 0, safety_min_n %% 1 == 0)
  list(n_dose = as.integer(n_dose), target = target, ei = ei, xi = xi,
       safety_min_n = as.integer(safety_min_n))
}

# Appendix pp. 16-17: p | x,n ~ Beta(1+x, 1+n-x).
overdose_probability <- function(x, n, target = 0.30) {
  stopifnot(all(is.finite(x)), all(is.finite(n)),
            all(x >= 0), all(n >= x), all(x %% 1 == 0), all(n %% 1 == 0))
  pbeta(target, 1 + x, 1 + n - x, lower.tail = FALSE)
}

# Split both tails into intervals of the EI width. Only the intervals
# touching 0 or 1 can be shorter; divide their mass by their ACTUAL width.
mtpi2_intervals <- function(ei = c(0.25, 0.35)) {
  width <- diff(ei)
  lower <- ei[1] - (0:ceiling(ei[1] / width)) * width
  upper <- ei[2] + (0:ceiling((1 - ei[2]) / width)) * width
  cuts <- sort(unique(round(c(0, lower[lower > 0],
                              upper[upper < 1], 1), 12)))
  lo <- head(cuts, -1)
  hi <- tail(cuts, -1)
  action <- ifelse(hi <= ei[1] + 1e-10, "E",
                   ifelse(lo >= ei[2] - 1e-10, "D", "S"))
  data.frame(lower = lo, upper = hi, action = action)
}

mtpi2_scores <- function(x, n, settings = basyc_settings()) {
  stopifnot(length(x) == 1, length(n) == 1)
  overdose_probability(x, n, settings$target) # Also validates the counts.
  tab <- mtpi2_intervals(settings$ei)
  tab$probability <- pbeta(tab$upper, 1 + x, 1 + n - x) -
                     pbeta(tab$lower, 1 + x, 1 + n - x)
  tab$upm <- tab$probability / (tab$upper - tab$lower)
  tab
}

# Return E=escalate, S=stay, D=de-escalate, DU=de-escalate and exclude.
# The appendix gives NO minimum n for its posterior safety rule (default 0).
# Set safety_min_n=3 ONLY to reproduce Figure 2's small-sample convention.
# The current simulation passes 4: strict n > 3, a separate hypothesis.
mtpi2_decision <- function(x, n, settings = basyc_settings()) {
  risk <- overdose_probability(x, n, settings$target)
  if (n >= settings$safety_min_n && risk > settings$xi) return("DU")
  if (n == 0) return("S") # Explicit convention; no informed movement yet.
  tab <- mtpi2_scores(x, n, settings)
  winners <- tab$action[abs(tab$upm - max(tab$upm)) < 1e-12]
  # Exact UPM ties: prefer the more conservative action. No such tie
  # affects the Table 2 replay; this is an explicit numerical convention.
  c("D", "S", "E")[c("D", "S", "E") %in% winners][1]
}

decision_tables <- function(n_max = 15L, settings = basyc_settings()) {
  rows <- lapply(seq_len(n_max), function(n) {
    x <- 0:n
    risk <- overdose_probability(x, n, settings$target)
    data.frame(n = n, x = x, overdose = risk,
               DCR = ifelse(n >= settings$safety_min_n &
                              risk > settings$xi, "NG", "G"),
               DAR = vapply(x, mtpi2_decision, character(1),
                            n = n, settings = settings))
  })
  do.call(rbind, rows)
}
