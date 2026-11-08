# Data and presentation helpers. These do not participate in dose decisions.

# Run from the shared project root or from the Lyu_2018 folder.
# All exercise inputs and outputs stay inside Lyu_2018 in either case.
lyu_path <- function(...) {
  folder <- if (file.exists("Lyu_2018/R/mtpi2.R")) "Lyu_2018" else "."
  file.path(folder, ...)
}

read_table2_outcomes <- function() {
  d <- read.csv(lyu_path("data", "table2_outcomes.csv"))
  as.matrix(d[, c("y1", "y2", "y3")])
}

compare_table2_doses <- function(fit) {
  expected <- read.csv(lyu_path("data", "table2_expected_doses.csv"))
  a <- fit$assignments
  a$expected <- as.matrix(expected[, c("d1", "d2", "d3")])[
    cbind(match(a$patient, expected$patient), a$cycle)]
  a$match <- a$dose == a$expected
  a[order(a$patient, a$cycle), c("patient", "cohort", "cycle", "week",
                               "dose", "expected", "match")]
}

compare_table2_totals <- function(fit) {
  printed <- read.csv(lyu_path("data", "table2_printed_totals.csv"), check.names = FALSE)
  rows <- list()
  for (r in seq_len(nrow(printed))) {
    week <- printed$week[r]
    snap <- fit$snapshots[[as.character(week)]]
    for (j in seq_len(fit$n_cycle)) {
      for (d in seq_len(fit$settings$n_dose)) {
        name <- paste0("C", j, "_D", d)
        actual <- paste0(snap$x[d, j], "/", snap$n[d, j])
        rows[[length(rows) + 1]] <- data.frame(week, cycle = j, dose = d,
                         printed = printed[[name]][r], computed = actual,
                         match = actual == printed[[name]][r])
      }
    }
  }
  do.call(rbind, rows)
}

# A compact, patient-level equivalent of a Table 2 panel.
# O/X is an OBSERVED outcome; 'Assigned D#' is a FUTURE cycle assignment.
table2_at_week <- function(fit, week) {
  a <- fit$assignments[fit$assignments$week <= week, ]
  ids <- sort(unique(a$patient))
  out <- data.frame(Cohort = (ids - 1L) %/% fit$cohort_size + 1L, Patient = ids)
  for (j in seq_len(fit$n_cycle)) {
    out[[paste0("Cycle ", j)]] <- vapply(ids, function(id) {
      r <- a[a$patient == id & a$cycle == j, ]
      if (!nrow(r)) return("")
      if (r$end_week > week) return(paste("Assigned D", r$dose, sep = ""))
      y <- fit$observed$tox[fit$observed$patient == id & fit$observed$cycle == j]
      paste0(if (y == 1) "X" else "O", " at D", r$dose)
    }, character(1))
  }
  out$Status <- vapply(ids, function(id) {
    r <- a[a$patient == id, ]
    if (any(r$cycle == fit$n_cycle & r$end_week <= week)) return("Exit")
    pending <- r[r$end_week > week, ]
    if (nrow(pending)) {
      if (pending$cycle == 1) "Enrolled" else paste("Go to cycle", pending$cycle)
    } else "Stopped"
  }, character(1))
  out
}

# Independent manual transcription of Figure 2B, columns n=1,...,15.
# Each string lists actions for x=0,...,n. The n=2,x=2 entry is D,
# although the appendix's unqualified posterior rule would give DU.
figure2_reference <- function() {
  columns <- c(
    "E D", "E D D", "E S D DU", "E S D DU DU",
    "E E D D DU DU", "E E S D DU DU DU", "E E S D D DU DU DU",
    "E E S D D DU DU DU DU", "E E E S D DU DU DU DU DU",
    "E E E S D D DU DU DU DU DU", "E E E S D D DU DU DU DU DU DU",
    "E E E S S D D DU DU DU DU DU DU",
    "E E E E S D D DU DU DU DU DU DU DU",
    "E E E E S D D D DU DU DU DU DU DU DU",
    "E E E E S S D D DU DU DU DU DU DU DU DU")
  # Figure 2A: first x at which the decision is NG; Inf = no NG in column.
  first_ng <- c(Inf, Inf, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8)
  do.call(rbind, lapply(seq_along(columns), function(n) {
    actions <- strsplit(columns[n], " ", fixed = TRUE)[[1]]
    stopifnot(length(actions) == n + 1)
    data.frame(n = n, x = 0:n, DAR = actions,
               DCR = ifelse(0:n >= first_ng[n], "NG", "G"))
  }))
}

# Plot the two small decision modules in the format of Figure 2.
plot_decision_tables <- function(tables) {
  old <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  colors <- c(G = "#c8dfec", NG = "#dfb4ac", E = "#c8dfec",
               S = "#f4e2a7", D = "#c1b57c", DU = "#dfb4ac")
  max_n <- max(tables$n)
  for (kind in c("DCR", "DAR")) {
    plot.new()
    plot.window(xlim = c(0.5, max_n + 0.5), ylim = c(max_n + 0.5, -0.5))
    for (r in seq_len(nrow(tables))) {
      n <- tables$n[r]
      x <- tables$x[r]
      label <- tables[[kind]][r]
      rect(n - 0.5, x - 0.5, n + 0.5, x + 0.5, col = colors[label], border = "white")
      text(n, x, label, cex = 0.65)
    }
    axis(1, at = seq_len(max_n))
    axis(2, at = 0:max_n, las = 1)
    title(main = kind, xlab = "Number of observed patients (n)", ylab = "Number of DLTs (x)")
  }
}
