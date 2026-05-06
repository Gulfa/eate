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
  source("../../det_model.R", chdir = TRUE)
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
# get_eate_network
# ---------------------------------------------------------------------------

# Small fixed contact matrix reused across all get_eate_network tests
local({
  set.seed(99)
  N_net <- 20; t_net <- 5; n_vac_net <- 3
  c_ij_net <- matrix(0, N_net, N_net)
  for (i in 1:(N_net-1)) for (j in (i+1):N_net)
    c_ij_net[i,j] <- c_ij_net[j,i] <- rbinom(1, 1, 6/N_net)
  for (i in 1:N_net) if (sum(c_ij_net[i,]) == 0) {
    j <- sample(setdiff(1:N_net, i), 1); c_ij_net[i,j] <- c_ij_net[j,i] <- 1
  }
  k_mean_net <- mean(rowSums(c_ij_net))
  assign("N_net",    N_net,    envir = parent.env(environment()))
  assign("t_net",    t_net,    envir = parent.env(environment()))
  assign("n_vac_net",n_vac_net,envir = parent.env(environment()))
  assign("c_ij_net", c_ij_net, envir = parent.env(environment()))
  assign("k_mean_net",k_mean_net,envir=parent.env(environment()))
})

test_that("get_eate_network method='full' returns correct structure", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")
  set.seed(1)
  res <- get_eate_network(alpha=0.5, beta=1, f=0.5, N=N_net, t=t_net,
                          c_ij=c_ij_net, n_vac=n_vac_net, method="full",
                          k_mean=k_mean_net, mc.cores=1)
  expect_named(res, c("t", "eate", "num", "denom", "method", "sim"))
  expect_setequal(unique(res$method), c("full", "CRR"))
  expect_equal(nrow(res), n_vac_net * t_net * 2)
  expect_setequal(unique(res$t), 1:t_net)
})

test_that("get_eate_network method='frozen' returns correct structure", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")
  skip_if_not_installed("deSolve")
  set.seed(2)
  res <- get_eate_network(alpha=0.5, beta=1, f=0.5, N=N_net, t=t_net,
                          c_ij=c_ij_net, n_vac=n_vac_net, method="frozen",
                          k_mean=k_mean_net, slowdown=2, mc.cores=1)
  expect_named(res, c("t", "eate", "num", "denom", "method", "sim"))
  expect_setequal(unique(res$method), c("frozen", "CRR"))
  expect_equal(nrow(res), n_vac_net * t_net * 2)
  expect_setequal(unique(res$t), 1:t_net)
})

test_that("get_eate_network method='both' returns correct structure", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")
  skip_if_not_installed("deSolve")
  set.seed(3)
  res <- get_eate_network(alpha=0.5, beta=1, f=0.5, N=N_net, t=t_net,
                          c_ij=c_ij_net, n_vac=n_vac_net, method="both",
                          k_mean=k_mean_net, slowdown=2, mc.cores=1)
  expect_named(res, c("t", "eate", "num", "denom", "method", "sim"))
  expect_setequal(unique(res$method), c("full", "frozen", "CRR"))
  expect_equal(nrow(res), n_vac_net * t_net * 3)
  expect_setequal(unique(res$t), 1:t_net)
  # one CRR row per sim (not duplicated between full and frozen)
  expect_equal(nrow(res[res$method == "CRR", ]), n_vac_net * t_net)
})

test_that("get_eate_network eate and CRR values are finite and positive", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")
  set.seed(4)
  res <- get_eate_network(alpha=0.5, beta=1, f=0.5, N=N_net, t=t_net,
                          c_ij=c_ij_net, n_vac=n_vac_net, method="full",
                          k_mean=k_mean_net, mc.cores=1)
  expect_true(all(is.finite(res$eate)))
  expect_true(all(res$eate[res$method == "CRR"] > 0))
})

# ---------------------------------------------------------------------------
# get_frailty_eate
# ---------------------------------------------------------------------------

test_that("get_frailty_eate method='full' returns correct structure", {
  set.seed(10)
  res <- get_frailty_eate(alpha=0.5, sd=0.2, beta=1, f=0.5, N=100, t=5,
                          n_frailty=5, method="full", n_vac=3, mc.cores=1)
  expect_named(res, c("t", "eate", "num", "denom", "method", "sim"))
  expect_setequal(unique(res$method), c("full", "CRR"))
  expect_equal(nrow(res), 3 * 5 * 2)
  expect_setequal(unique(res$t), 1:5)
})

test_that("get_frailty_eate method='frozen' returns correct structure", {
  set.seed(11)
  res <- get_frailty_eate(alpha=0.5, sd=0.2, beta=1, f=0.5, N=100, t=5,
                          n_frailty=5, method="frozen", slowdown=2, n_vac=3, mc.cores=1)
  expect_named(res, c("t", "eate", "num", "denom", "method", "sim"))
  expect_setequal(unique(res$method), c("frozen", "CRR"))
  expect_equal(nrow(res), 3 * 5 * 2)
  expect_setequal(unique(res$t), 1:5)
})

test_that("get_frailty_eate method='both' returns correct structure", {
  set.seed(12)
  res <- get_frailty_eate(alpha=0.5, sd=0.2, beta=1, f=0.5, N=100, t=5,
                          n_frailty=5, method="both", slowdown=2, n_vac=3, mc.cores=1)
  expect_named(res, c("t", "eate", "num", "denom", "method", "sim"))
  expect_setequal(unique(res$method), c("full", "frozen", "CRR"))
  expect_equal(nrow(res), 3 * 5 * 3)
  expect_setequal(unique(res$t), 1:5)
  # one CRR set per replication (not duplicated)
  expect_equal(nrow(res[res$method == "CRR", ]), 3 * 5)
})

test_that("get_frailty_eate eate and CRR values are finite and positive at late times", {
  set.seed(13)
  res <- get_frailty_eate(alpha=0.5, sd=0.2, beta=1, f=0.5, N=500, t=20,
                          n_frailty=10, method="full", n_vac=3, mc.cores=1)
  # once the epidemic is established both EATE and CRR should be finite and positive
  late <- res[res$t >= 10, ]
  expect_true(all(is.finite(late$eate)))
  expect_true(all(late$eate[late$method == "CRR"] > 0))
})

# ---------------------------------------------------------------------------
# run_det_ncd
# ---------------------------------------------------------------------------

test_that("run_det_ncd returns full_results with correct dimensions", {
  res <- run_det_ncd(N = c(500, 500), t = 20, susceptibility = c(1, 0.5))
  expect_true(!is.null(res$full_results))
  expect_equal(nrow(res$full_results), 20)
})
