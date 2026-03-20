# Tests for model.R
# Run from project root so odin can find det_mod_cd.R / det_mod_ncd.R
#
# Requires: odin, dplyr, data.table

skip_if_not_installed("odin")
skip_if_not_installed("dplyr")
skip_if_not_installed("data.table")
skip_if_not_installed("Pareto")

library(data.table)
library(dplyr)

# Source once per session; guard against re-sourcing
if (!exists("det_model_cd")) {
  # When run via test_dir(), working dir is tests/testthat; ../../ is project root
  source("../../model.R", chdir = TRUE)
}

# ---------------------------------------------------------------------------
# get_frailty
# ---------------------------------------------------------------------------

test_that("get_frailty returns a valid discrete probability distribution", {
  f <- get_frailty(mean = 0.5, sd = 0.2, n = 100)
  expect_length(f$x, 100)
  expect_length(f$p, 100)
  expect_equal(sum(f$p), 1, tolerance = 1e-10)
  expect_true(all(f$p >= 0))
  expect_true(all(f$x > 0) && all(f$x < 1))
})

test_that("get_frailty x values are ordered and within (0,1)", {
  f <- get_frailty(mean = 0.3, sd = 0.1, n = 50)
  expect_true(all(diff(f$x) > 0))          # strictly increasing
  expect_true(min(f$x) > 0)
  expect_true(max(f$x) < 1)
})

test_that("get_frailty with larger sd gives more spread", {
  f_narrow <- get_frailty(mean = 0.5, sd = 0.05, n = 200)
  f_wide   <- get_frailty(mean = 0.5, sd = 0.30, n = 200)
  # Weighted variance should be larger for wider distribution
  wvar <- function(f) sum(f$p * (f$x - 0.5)^2)
  expect_gt(wvar(f_wide), wvar(f_narrow))
})

# ---------------------------------------------------------------------------
# cij_NGM
# ---------------------------------------------------------------------------

test_that("cij_NGM returns positive beta_R", {
  N   <- c(500, 500)
  cij <- matrix(c(1, 1, 1, 1), nrow = 2) / 2
  res <- cij_NGM(cij, N, susceptibility = c(1, 1),
                 transmisibility = c(1, 1), gamma = 1/3)
  expect_true(is.numeric(res$beta_R))
  expect_gt(res$beta_R, 0)
})

test_that("cij_NGM beta_R scales with susceptibility", {
  N   <- c(500, 500)
  cij <- matrix(c(1, 1, 1, 1), nrow = 2) / 2
  res_full <- cij_NGM(cij, N, susceptibility = c(1, 1),
                      transmisibility = c(1, 1), gamma = 1/3)
  res_half <- cij_NGM(cij, N, susceptibility = c(0.5, 0.5),
                      transmisibility = c(1, 1), gamma = 1/3)
  expect_equal(res_half$beta_R, res_full$beta_R * 0.5, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# get_beta
# ---------------------------------------------------------------------------

test_that("get_beta returns a positive scalar", {
  beta <- get_beta(R = 2, alpha = 0.5, sd = 0.2, f = 0.5,
                   N = 500, n_frailty = 50, gamma = 1/2)
  expect_length(beta, 1)
  expect_gt(beta, 0)
})

test_that("get_beta is monotone in R", {
  beta2 <- get_beta(R = 2, alpha = 0.5, sd = 0.2, f = 0.5, N = 500, n_frailty = 50)
  beta3 <- get_beta(R = 3, alpha = 0.5, sd = 0.2, f = 0.5, N = 500, n_frailty = 50)
  expect_gt(beta3, beta2)
})

# ---------------------------------------------------------------------------
# run_det_cd
# ---------------------------------------------------------------------------

test_that("run_det_cd returns expected structure", {
  mm    <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)
  N     <- c(500, 500)
  t     <- 30
  I_ini <- c(2, 0)
  res <- run_det_cd(mm, rep(0.5, t), N, t, I_ini,
                    susceptibility = c(1, 0.5), gamma = 1/3)
  expect_named(res, c("main", "ind", "full_results"))
  expect_equal(nrow(res$main), t)
  expect_true(all(c("t", "CRR", "HR", "unvac", "vac") %in% colnames(res$main)))
})

test_that("run_det_cd cumulative cases are non-decreasing", {
  mm    <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)
  N     <- c(500, 500)
  t     <- 50
  I_ini <- c(5, 0)
  res <- run_det_cd(mm, rep(0.5, t), N, t, I_ini,
                    susceptibility = c(1, 0.5), gamma = 1/3)
  expect_true(all(diff(res$main$unvac) >= -1e-6))
  expect_true(all(diff(res$main$vac)   >= -1e-6))
})

test_that("run_det_cd CRR ~ 1 when both groups have identical susceptibility", {
  mm    <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)
  N     <- c(500, 500)
  t     <- 80
  I_ini <- c(5, 5)
  res <- run_det_cd(mm, rep(0.5, t), N, t, I_ini,
                    susceptibility = c(1, 1), gamma = 1/3)
  final_CRR <- tail(res$main$CRR, 1)
  expect_equal(final_CRR, 1, tolerance = 0.02)
})

test_that("run_det_cd vaccination (alpha < 1) reduces cumulative attack rate", {
  mm    <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)
  N     <- c(1000, 1000)
  t     <- 100
  I_ini <- c(10, 0)
  res <- run_det_cd(mm, rep(0.5, t), N, t, I_ini,
                    susceptibility = c(1, 0.3), gamma = 1/3)
  final_CRR <- tail(res$main$CRR, 1)
  expect_lt(final_CRR, 1)
})

# ---------------------------------------------------------------------------
# sparse vs dense equivalence
# ---------------------------------------------------------------------------

test_that("sparse run_det_cd matches dense on a binary contact matrix", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  # 4-node binary network: nodes 1,2 unvaccinated; 3,4 vaccinated
  cm <- matrix(c(0,1,1,0,
                 1,0,0,1,
                 1,0,0,1,
                 0,1,1,0), nrow=4, byrow=TRUE)
  N       <- c(500, 500, 500, 500)
  t       <- 50
  beta_day <- rep(0.5, t)
  I_ini   <- c(2, 0, 0, 0)
  sus     <- c(1, 1, 0.5, 0.5)

  res_dense  <- run_det_cd(cm, beta_day, N, t, I_ini,
                           susceptibility = sus, gamma = 1/3, sparse = FALSE)
  res_sparse <- run_det_cd(cm, beta_day, N, t, I_ini,
                           susceptibility = sus, gamma = 1/3, sparse = TRUE)

  # CRR should agree closely (ODE solvers differ, so allow small tolerance)
  expect_equal(tail(res_sparse$main$CRR, 1),
               tail(res_dense$main$CRR,  1),
               tolerance = 1e-3)

  # Cumulative cases should agree within 1 person
  expect_equal(res_sparse$main$vac,   res_dense$main$vac,   tolerance = 1)
  expect_equal(res_sparse$main$unvac, res_dense$main$unvac, tolerance = 1)
})

# ---------------------------------------------------------------------------
# run_det_ncd
# ---------------------------------------------------------------------------

test_that("run_det_ncd returns full_results with correct dimensions", {
  res <- run_det_ncd(N = c(500, 500), t = 20, susceptibility = c(1, 0.5))
  expect_true(!is.null(res$full_results))
  expect_equal(nrow(res$full_results), 20)
})
