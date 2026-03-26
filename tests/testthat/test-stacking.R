# Tests for stacking.R
# These tests do NOT require odin and can run independently.
#
# Requires: dplyr, data.table, zoo, adaptivetau, parallel

skip_if_not_installed("dplyr")
skip_if_not_installed("data.table")
skip_if_not_installed("zoo")
skip_if_not_installed("adaptivetau")

if (!exists("run_stacking")) {
  # When run via test_dir(), working dir is tests/testthat; ../../ is project root
  source("../../stoch_model.R")
}

library(dplyr)
library(data.table)

# ---------------------------------------------------------------------------
# regularise
# ---------------------------------------------------------------------------

test_that("regularise returns exactly the requested timepoints", {
  # Input times must not overlap with timepoints to avoid duplicate rows
  df <- data.frame(N = c(10, 12, 15), time = c(0.5, 1.5, 2.5))
  timepoints <- c(0, 1, 2, 3)
  result <- regularise(df, timepoints)
  expect_setequal(result$time, timepoints)
  expect_equal(nrow(result), length(timepoints))
})

test_that("regularise carries last observation forward (LOCF)", {
  # Input times must not overlap with timepoints
  df <- data.frame(N = c(10, 20), time = c(0.5, 1.5))
  result <- regularise(df, c(0, 1, 2))
  # At t=1 (after obs at 0.5), LOCF gives 10
  expect_equal(result$N[result$time == 1], 10)
  # At t=2 (after obs at 1.5), LOCF gives 20
  expect_equal(result$N[result$time == 2], 20)
})

# ---------------------------------------------------------------------------
# generate_contacts
# ---------------------------------------------------------------------------

test_that("generate_contacts returns a list with events before T", {
  set.seed(1)
  contacts <- generate_contacts(N = 8, T = 2, contact_rate = 0.5)
  expect_true(is.list(contacts))
  if (length(contacts) > 0) {
    tvals <- sapply(contacts, function(x) x$t)
    expect_true(all(tvals < 2))
    types <- sapply(contacts, function(x) x$type)
    expect_true(all(types == "contact"))
    froms <- sapply(contacts, function(x) x$from)
    tos   <- sapply(contacts, function(x) x$to)
    # No self-contacts
    expect_true(all(froms != tos))
    expect_true(all(froms >= 1 & froms <= 8))
    expect_true(all(tos   >= 1 & tos   <= 8))
  }
})

# ---------------------------------------------------------------------------
# generate_linear_event_times
# ---------------------------------------------------------------------------

test_that("generate_linear_event_times returns events with t < T", {
  set.seed(2)
  events <- generate_linear_event_times(N = 20, rate = 0.5, T = 3)
  expect_true(is.list(events))
  if (length(events) > 0) {
    tvals <- sapply(events, function(x) x$t)
    expect_true(all(tvals < 3))
    whos  <- sapply(events, function(x) x$who)
    unifs <- sapply(events, function(x) x$unif)
    expect_true(all(whos >= 1 & whos <= 20))
    expect_true(all(unifs >= 0 & unifs <= 1))
  }
})

test_that("generate_linear_event_times average count matches Poisson expectation", {
  set.seed(3)
  N <- 50; rate <- 0.2; T <- 5
  # Expected events per individual ~ rate * T, total ~ N * rate * T
  counts <- replicate(200, length(generate_linear_event_times(N, rate, T)))
  expect_equal(mean(counts), N * rate * T, tolerance = N * rate * T * 0.15)
})

# ---------------------------------------------------------------------------
# run_events_linear
# ---------------------------------------------------------------------------

test_that("run_events_linear: no events leaves state unchanged", {
  state <- rep(0, 10)
  result <- run_events_linear(state, rep(1, 10), list(), T = 5)
  expect_equal(result$last_state, state)
})

test_that("run_events_linear: state only increases (0 -> 1 transitions)", {
  set.seed(4)
  state  <- rep(0, 20)
  sus    <- rep(1, 20)
  events <- generate_linear_event_times(20, rate = 0.5, T = 4)
  result <- run_events_linear(state, sus, events, T = 4)
  # State is binary: 0 (susceptible) or 1 (infected)
  expect_true(all(result$last_state %in% c(0, 1)))
  # Infections are non-negative
  expect_gte(sum(result$last_state), 0)
})

test_that("run_events_linear: zero susceptibility prevents infection", {
  set.seed(5)
  state  <- rep(0, 10)
  sus    <- rep(0, 10)   # alpha = 0 for everyone
  events <- generate_linear_event_times(10, rate = 1, T = 5)
  result <- run_events_linear(state, sus, events, T = 5)
  expect_equal(sum(result$last_state), 0)
})

# ---------------------------------------------------------------------------
# run_stacking
# ---------------------------------------------------------------------------

# Minimal birth-death rate function used in the original test() function
bd_rate <- function(x, params, t) {
  c(params$beta * x["N"], params$delta * x["N"])
}

test_that("run_stacking returns X and Z with identical time vectors", {
  set.seed(6)
  init <- c(N = 50)
  trans <- list(c(N = 1), c(N = -1))
  p <- list(beta = 0.03, delta = 0.02)
  result <- run_stacking(init, init, trans, p, p, bd_rate, bd_rate, T = 5)
  expect_named(result, c("X", "Z"))
  expect_equal(result$X$time, result$Z$time)
})

test_that("run_stacking with identical params produces identical X and Z", {
  set.seed(7)
  init <- c(N = 30)
  trans <- list(c(N = 1), c(N = -1))
  p <- list(beta = 0.04, delta = 0.03)
  result <- run_stacking(init, init, trans, p, p, bd_rate, bd_rate, T = 3)
  # With same params and shared coupling, X and Z must be identical
  expect_equal(result$X$N, result$Z$N)
})

test_that("run_stacking terminates and first row matches initial values", {
  set.seed(8)
  init1 <- c(N = 40)
  init2 <- c(N = 60)
  trans <- list(c(N = 1), c(N = -1))
  p1 <- list(beta = 0.03, delta = 0.02)
  p2 <- list(beta = 0.05, delta = 0.02)
  result <- run_stacking(init1, init2, trans, p1, p2, bd_rate, bd_rate, T = 4)
  expect_equal(result$X$N[1], 40)
  expect_equal(result$Z$N[1], 60)
})

# ---------------------------------------------------------------------------
# Informal validation: linear implementation vs adaptivetau
# (based on test_linear_implementation in stacking.R)
# Checks that both approaches produce similar means under many simulations.
# ---------------------------------------------------------------------------

test_that("custom linear simulator mean matches adaptivetau over many sims", {
  skip_if_not_installed("adaptivetau")
  set.seed(9)

  N_sim <- 500
  N     <- 40
  T     <- 3
  rate  <- 0.08
  sus   <- c(rep(1, N / 2), rep(0.5, N / 2))
  state <- rep(0, N)

  custom_means <- unlist(parallel::mclapply(1:N_sim, function(i) {
    evts <- generate_linear_event_times(N, rate = rate, T = T)
    sum(run_events_linear(state, sus, evts, T)$last_state[1:(N / 2)])
  }, mc.cores = 1))

  LinearRateF <- function(x, p, t) {
    c(p$beta * x["S1"], 0, p$alpha * p$beta * x["S2"], 0)
  }
  init.vals <- c(S1 = N / 2, I1 = 0, R1 = 0, C1 = 0,
                 S2 = N / 2, I2 = 0, R2 = 0, C2 = 0)
  transitions <- list(c(S1 = -1, I1 = 1, C1 = 1), c(I1 = -1, R1 = 1),
                      c(S2 = -1, I2 = 1, C2 = 1), c(I2 = -1, R2 = 1))
  ssa_means <- unlist(parallel::mclapply(1:N_sim, function(i) {
    a <- adaptivetau::ssa.adaptivetau(init.vals, transitions, LinearRateF,
                                      list(beta = rate, alpha = 0.5), tf = T)
    a[nrow(a), "C1"]
  }, mc.cores = 1))

  # Means should agree within 15% of the larger mean
  tol <- max(mean(custom_means), mean(ssa_means)) * 0.15
  expect_equal(mean(custom_means), mean(ssa_means), tolerance = tol)
})
