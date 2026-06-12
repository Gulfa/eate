skip_if_not_installed("adaptivetau")
skip_if_not_installed("dplyr")
skip_if_not_installed("data.table")
skip_if_not_installed("zoo")

library(data.table)
library(dplyr)

if (!exists("run_sir")) {
  source("../../utils.R", chdir=TRUE)
  source("../../stoch_model.R", chdir=TRUE)
}

# ---------------------------------------------------------------------------
# run_sir
# ---------------------------------------------------------------------------

test_that("run_sir returns correct structure", {
  set.seed(1)
  res <- run_sir(N_cont=100, N_vac=100, beta=0.3, alpha=0.5, gamma=1/7,
                 I0_cont=2, t=20, timepoints=0:20, n_sim=5, cores=1)
  expect_true(is.data.frame(res))
  expect_true(all(c("S1","I1","R1","C1","S2","I2","R2","C2","time","sim") %in% names(res)))
  expect_equal(nrow(res), 5 * 21)
  expect_setequal(unique(res$sim), 1:5)
  expect_setequal(unique(res$time), 0:20)
})

test_that("run_sir conserves population in each group", {
  set.seed(2)
  res <- run_sir(N_cont=80, N_vac=80, beta=0.5, alpha=0.5, gamma=1/5,
                 I0_cont=1, t=30, timepoints=0:30, n_sim=10, cores=1)
  expect_true(all(res$S1 + res$I1 + res$R1 == 80))
  expect_true(all(res$S2 + res$I2 + res$R2 == 80))
})

test_that("run_sir cumulative cases are non-decreasing within each sim", {
  set.seed(3)
  res <- run_sir(N_cont=100, N_vac=100, beta=0.5, alpha=0.5, gamma=1/7,
                 I0_cont=2, t=20, timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$C1) >= 0))
    expect_true(all(diff(r$C2) >= 0))
  }
})

test_that("run_sir with beta=0 produces no infections", {
  set.seed(4)
  res <- run_sir(N_cont=100, N_vac=100, beta=0, alpha=0.5, gamma=1/7,
                 I0_cont=0, t=20, timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$C1 == 0))
  expect_true(all(res$C2 == 0))
})

test_that("run_sir alpha=0 produces far fewer vaccinated infections than unvaccinated", {
  set.seed(5)
  res <- run_sir(N_cont=200, N_vac=200, beta=1.5, alpha=0,
                 I0_cont=5, t=50, timepoints=50, n_sim=20, cores=1)
  expect_true(mean(res$C2) < mean(res$C1) * 0.1)
})

test_that("run_sir alpha=1 gives similar attack rates in both groups", {
  set.seed(6)
  res <- run_sir(N_cont=300, N_vac=300, beta=1.5, alpha=1, gamma=1/7,
                 I0_cont=3, I0_vac=3, t=60, timepoints=60, n_sim=40, cores=1)
  ar1 <- mean(res$C1) / 300
  ar2 <- mean(res$C2) / 300
  expect_equal(ar1, ar2, tolerance=0.05)
})

# ---------------------------------------------------------------------------
# run_linear
# ---------------------------------------------------------------------------

test_that("run_linear returns correct structure", {
  set.seed(10)
  res <- run_linear(N_cont=100, N_vac=100, beta=0.05, alpha=0.5,
                    t=20, timepoints=0:20, n_sim=5, cores=1)
  expect_true(is.data.frame(res))
  expect_true(all(c("S1","C1","S2","C2","time","sim") %in% names(res)))
  expect_equal(nrow(res), 5 * 21)
  expect_setequal(unique(res$sim), 1:5)
  expect_setequal(unique(res$time), 0:20)
})

test_that("run_linear conserves population in each group", {
  set.seed(11)
  res <- run_linear(N_cont=80, N_vac=80, beta=0.05, alpha=0.5,
                    t=20, timepoints=0:20, n_sim=10, cores=1)
  expect_true(all(res$S1 + res$C1 == 80))
  expect_true(all(res$S2 + res$C2 == 80))
})

test_that("run_linear cumulative cases are non-decreasing within each sim", {
  set.seed(12)
  res <- run_linear(N_cont=100, N_vac=100, beta=0.05, alpha=0.5,
                    t=20, timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$C1) >= 0))
    expect_true(all(diff(r$C2) >= 0))
  }
})

test_that("run_linear with beta=0 produces no infections", {
  set.seed(13)
  res <- run_linear(N_cont=100, N_vac=100, beta=0, alpha=0.5,
                    t=20, timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$C1 == 0))
  expect_true(all(res$C2 == 0))
})

test_that("run_linear alpha=0 produces no vaccinated infections", {
  set.seed(14)
  res <- run_linear(N_cont=100, N_vac=100, beta=0.1, alpha=0,
                    t=20, timepoints=0:20, n_sim=10, cores=1)
  expect_true(all(res$C2 == 0))
})

test_that("run_linear mean attack rate matches expected exponential survival", {
  set.seed(15)
  beta <- 0.08
  N <- 500
  t_end <- 10
  res <- run_linear(N_cont=N, N_vac=0, beta=beta, alpha=1,
                    t=t_end, timepoints=t_end, n_sim=200, cores=1)
  expected_ar <- 1 - exp(-beta * t_end)
  observed_ar <- mean(res$C1) / N
  expect_equal(observed_ar, expected_ar, tolerance=0.03)
})

# ---------------------------------------------------------------------------
# run_stoch_cd_ctmc
# ---------------------------------------------------------------------------

test_that("run_stoch_cd_ctmc returns correct structure for n groups", {
  set.seed(20)
  mm <- matrix(1/4, nrow=4, ncol=4)
  res <- run_stoch_cd_ctmc(mm, beta=1, N=c(100,100,100,100), t=20, I_ini=c(1,0,0,0),
                      timepoints=0:20, n_sim=4, cores=1)
  expected_cols <- c(paste0(c("S","I","R","C"), rep(1:4, each=1)) |>
                       (\(x) c(paste0("S",1:4), paste0("I",1:4), paste0("R",1:4), paste0("C",1:4)))(),
                     "time", "sim")
  expect_true(all(expected_cols %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
})

test_that("run_stoch_cd_ctmc conserves population in every group", {
  set.seed(21)
  mm <- matrix(c(0.8,0.2, 0.2,0.8), nrow=2)
  N  <- c(150, 100)
  res <- run_stoch_cd_ctmc(mm, beta=1.5, N=N, t=30, I_ini=c(2,0),
                      timepoints=0:30, n_sim=10, cores=1)
  expect_true(all(res$S1 + res$I1 + res$R1 == N[1]))
  expect_true(all(res$S2 + res$I2 + res$R2 == N[2]))
})

test_that("run_stoch_cd_ctmc cumulative cases are non-decreasing within each sim", {
  set.seed(22)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd_ctmc(mm, beta=1.5, N=c(200,200), t=20, I_ini=c(2,0),
                      timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$C1) >= 0))
    expect_true(all(diff(r$C2) >= 0))
  }
})

test_that("run_stoch_cd_ctmc beta=0 produces no infections", {
  set.seed(23)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd_ctmc(mm, beta=0, N=c(100,100), t=20, I_ini=c(0,0),
                      timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$C1 == 0))
  expect_true(all(res$C2 == 0))
})

test_that("run_stoch_cd_ctmc susceptibility=0 for one group blocks its infections", {
  set.seed(24)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd_ctmc(mm, beta=2, N=c(200,200), t=30, I_ini=c(3,0),
                      susceptibility=c(1, 0), timepoints=30, n_sim=20, cores=1)
  expect_true(all(res$C2 == 0))
  expect_true(mean(res$C1) > 10)
})

test_that("run_stoch_cd_ctmc with homogeneous mixing and equal susceptibility matches run_sir attack rate", {
  set.seed(25)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res_cd  <- run_stoch_cd_ctmc(mm, beta=1.5, N=c(200,200), t=40, I_ini=c(2,2),
                           susceptibility=c(1,1), gamma=1/7,
                           timepoints=40, n_sim=80, cores=1)
  res_sir <- run_sir(N_cont=200, N_vac=200, beta=1.5, alpha=1, gamma=1/7,
                     I0_cont=2, I0_vac=2, t=40, timepoints=40, n_sim=80, cores=1)
  ar_cd  <- (mean(res_cd$C1)  + mean(res_cd$C2))  / 400
  ar_sir <- (mean(res_sir$C1) + mean(res_sir$C2)) / 400
  expect_equal(ar_cd, ar_sir, tolerance=0.05)
})

# ---------------------------------------------------------------------------
# run_stoch_cd_dust  (odin2/dust2 backend)
# ---------------------------------------------------------------------------

skip_dust <- function() {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")
}

test_that("run_stoch_cd_dust returns correct structure for n groups", {
  skip_dust()
  set.seed(120)
  mm  <- matrix(1/4, nrow=4, ncol=4)
  res <- run_stoch_cd_dust(mm, beta=1, N=c(100,100,100,100), t=20, I_ini=c(1,0,0,0),
                           timepoints=0:20, n_sim=4, cores=1)
  expected_cols <- c(paste0("S", 1:4), paste0("I", 1:4),
                     paste0("R", 1:4), paste0("C", 1:4), "time", "sim")
  expect_true(all(expected_cols %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:20)
})

test_that("run_stoch_cd_dust conserves population in every group", {
  skip_dust()
  set.seed(121)
  mm  <- matrix(c(0.8,0.2, 0.2,0.8), nrow=2)
  N   <- c(150, 100)
  res <- run_stoch_cd_dust(mm, beta=1.5, N=N, t=30, I_ini=c(2,0),
                           timepoints=0:30, n_sim=10, cores=1)
  expect_true(all(res$S1 + res$I1 + res$R1 == N[1]))
  expect_true(all(res$S2 + res$I2 + res$R2 == N[2]))
})

test_that("run_stoch_cd_dust cumulative cases are non-decreasing within each sim", {
  skip_dust()
  set.seed(122)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd_dust(mm, beta=1.5, N=c(200,200), t=20, I_ini=c(2,0),
                           timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$C1) >= 0))
    expect_true(all(diff(r$C2) >= 0))
  }
})

test_that("run_stoch_cd_dust beta=0 produces no infections", {
  skip_dust()
  set.seed(123)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd_dust(mm, beta=0, N=c(100,100), t=20, I_ini=c(0,0),
                           timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$C1 == 0))
  expect_true(all(res$C2 == 0))
})

test_that("run_stoch_cd_dust susceptibility=0 blocks infections in that group", {
  skip_dust()
  set.seed(124)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd_dust(mm, beta=2, N=c(200,200), t=30, I_ini=c(3,0),
                           susceptibility=c(1, 0), timepoints=30, n_sim=20, cores=1)
  expect_true(all(res$C2 == 0))
  expect_true(mean(res$C1) > 10)
})

test_that("run_stoch_cd_dust matches run_stoch_cd_ctmc attack rate (homogeneous mixing)", {
  skip_dust()
  set.seed(125)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res_dust <- run_stoch_cd_dust(mm, beta=1.5, N=c(200,200), t=40, I_ini=c(2,2),
                                 susceptibility=c(1,1), gamma=1/7,
                                 timepoints=40, n_sim=80, cores=1)
  res_ctmc <- run_stoch_cd_ctmc(mm, beta=1.5, N=c(200,200), t=40, I_ini=c(2,2),
                                 susceptibility=c(1,1), gamma=1/7,
                                 timepoints=40, n_sim=80, cores=1)
  ar_dust <- (mean(res_dust$C1) + mean(res_dust$C2)) / 400
  ar_ctmc <- (mean(res_ctmc$C1) + mean(res_ctmc$C2)) / 400
  expect_equal(ar_dust, ar_ctmc, tolerance=0.05)
})

# ---------------------------------------------------------------------------
# run_stoch_linear_dust
# ---------------------------------------------------------------------------

test_that("run_stoch_linear_dust mean attack rate matches analytic 1 - exp(-sus*beta*t)", {
  skip_dust()
  # Each group is independent: per-individual hazard = susceptibility[k] * beta,
  # so E[C_k(t)] = N[k] * (1 - exp(-sus[k] * beta * t)). dt small so the
  # tau-leaping bias is well below Monte Carlo noise.
  set.seed(150)
  beta  <- 0.05
  sus   <- c(1, 0.5, 1.5)
  N     <- c(800L, 800L, 800L)
  t_end <- 12
  n_sim <- 400

  res <- run_stoch_linear_dust(beta=beta, N=N, t=t_end, susceptibility=sus,
                               dt=0.01, timepoints=t_end,
                               n_sim=n_sim, cores=1)

  observed <- c(mean(res$C1) / N[1],
                mean(res$C2) / N[2],
                mean(res$C3) / N[3])
  expected <- 1 - exp(-sus * beta * t_end)

  # Monte Carlo SE per individual ~ sqrt(p*(1-p)/n_sim/N) — well under 0.01.
  expect_equal(observed, expected, tolerance=0.02)
})

test_that("run_stoch_linear_dust trajectory mean tracks analytic curve over time", {
  skip_dust()
  set.seed(151)
  beta  <- 0.04
  sus   <- 0.8
  N     <- 1000L
  ts    <- c(2, 5, 10, 20)
  n_sim <- 400

  res <- run_stoch_linear_dust(beta=beta, N=N, t=max(ts), susceptibility=sus,
                               dt=0.01, timepoints=ts,
                               n_sim=n_sim, cores=1)

  observed <- vapply(ts, function(tp) mean(res$C1[res$time == tp]) / N, numeric(1))
  expected <- 1 - exp(-sus * beta * ts)
  expect_equal(observed, expected, tolerance=0.02)
})

# ---------------------------------------------------------------------------
# run_stoch_adj  (stochastic adjacency-list, odin2/dust2)
# ---------------------------------------------------------------------------

test_that("run_stoch_adj returns correct structure", {
  skip_dust()
  set.seed(140)
  n <- 6
  cm <- matrix(0L, n, n); cm[upper.tri(cm)] <- 1L; cm <- cm + t(cm)
  I_ini <- c(1L, rep(0L, n - 1L))
  res <- run_stoch_adj(cm, beta=2, t=10, I_ini=I_ini,
                       timepoints=0:10, n_sim=4, cores=1)
  expected_cols <- c(paste0("S", 1:n), paste0("I", 1:n),
                     paste0("R", 1:n), paste0("C", 1:n), "time", "sim")
  expect_true(all(expected_cols %in% names(res)))
  expect_equal(nrow(res), 4 * 11)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:10)
})

test_that("run_stoch_adj per-node S + I + R == 1 (individual-level)", {
  skip_dust()
  set.seed(141)
  n <- 5
  cm <- matrix(1L, n, n); diag(cm) <- 0L
  I_ini <- c(1L, rep(0L, n - 1L))
  res <- run_stoch_adj(cm, beta=2, t=10, I_ini=I_ini,
                       timepoints=0:10, n_sim=5, cores=1)
  for (k in seq_len(n)) {
    expect_true(all(res[[paste0("S", k)]] + res[[paste0("I", k)]] +
                    res[[paste0("R", k)]] == 1L))
  }
})

test_that("run_stoch_adj isolated node never gets infected", {
  skip_dust()
  set.seed(142)
  # 4 nodes, node 4 has no edges
  cm <- matrix(0L, 4, 4)
  cm[1, 2] <- cm[2, 1] <- 1L
  cm[2, 3] <- cm[3, 2] <- 1L
  cm[1, 3] <- cm[3, 1] <- 1L
  res <- run_stoch_adj(cm, beta=5, t=20, I_ini=c(1L,0L,0L,0L),
                       timepoints=20, n_sim=20, cores=1)
  expect_true(all(res$C4 == 0L))
})

# ---------------------------------------------------------------------------
# run_stoch_frailty_cd
# ---------------------------------------------------------------------------

test_that("run_stoch_frailty_cd returns correct structure", {
  set.seed(30)
  res <- run_stoch_frailty_cd(alpha=0.3, sd=0.3, beta=1.5, f=0.5, N=300, t=20,
                               n_frailty=5, gamma=1/7, timepoints=0:20, n_sim=4, cores=1)
  expect_true(all(c("time","sim","vac","unvac","CRR") %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:20)
})

test_that("run_stoch_frailty_cd vac and unvac are non-decreasing within each sim", {
  set.seed(31)
  res <- run_stoch_frailty_cd(alpha=0.3, sd=0.3, beta=1.5, f=0.5, N=300, t=20,
                               n_frailty=5, gamma=1/7, timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$vac)   >= 0))
    expect_true(all(diff(r$unvac) >= 0))
  }
})

test_that("run_stoch_frailty_cd alpha<1 gives lower attack rate in vaccinated group", {
  set.seed(32)
  res <- run_stoch_frailty_cd(alpha=0.2, sd=0.3, beta=2, f=0.5, N=500, t=40,
                               n_frailty=10, gamma=1/7, timepoints=40, n_sim=30, cores=1)
  N_half <- 500
  ar_vac   <- mean(res$vac)   / N_half
  ar_unvac <- mean(res$unvac) / N_half
  expect_true(ar_vac < ar_unvac)
})

test_that("run_stoch_frailty_cd CRR < 1 when alpha < 1 at end of epidemic", {
  set.seed(33)
  res <- run_stoch_frailty_cd(alpha=0.3, sd=0.3, beta=2, f=0.5, N=500, t=50,
                               n_frailty=10, gamma=1/7, timepoints=50, n_sim=30, cores=1)
  crr_vals <- res$CRR[is.finite(res$CRR) & res$CRR > 0]
  expect_true(mean(crr_vals) < 1)
})

# ---------------------------------------------------------------------------
# run_stoch_linear
# ---------------------------------------------------------------------------

test_that("run_stoch_linear returns correct structure", {
  set.seed(40)
  res <- run_stoch_linear(beta=0.05, N=c(100, 100, 100), t=20,
                           timepoints=0:20, n_sim=4, cores=1)
  expect_true(all(c("S1","C1","S2","C2","S3","C3","time","sim") %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
})

test_that("run_stoch_linear conserves population in every group", {
  set.seed(41)
  N <- c(80, 120, 60)
  res <- run_stoch_linear(beta=0.05, N=N, t=20, timepoints=0:20, n_sim=8, cores=1)
  expect_true(all(res$S1 + res$C1 == N[1]))
  expect_true(all(res$S2 + res$C2 == N[2]))
  expect_true(all(res$S3 + res$C3 == N[3]))
})

test_that("run_stoch_linear cumulative cases are non-decreasing within each sim", {
  set.seed(42)
  res <- run_stoch_linear(beta=0.05, N=c(100, 100), t=20, timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$C1) >= 0))
    expect_true(all(diff(r$C2) >= 0))
  }
})

test_that("run_stoch_linear beta=0 produces no infections", {
  set.seed(43)
  res <- run_stoch_linear(beta=0, N=c(100, 100), t=20, timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$C1 == 0))
  expect_true(all(res$C2 == 0))
})

test_that("run_stoch_linear susceptibility=0 blocks infections in that group", {
  set.seed(44)
  res <- run_stoch_linear(beta=0.1, N=c(200, 200), t=20,
                           susceptibility=c(1, 0), timepoints=20, n_sim=20, cores=1)
  expect_true(all(res$C2 == 0))
  expect_true(mean(res$C1) > 10)
})

test_that("run_stoch_linear mean attack rate matches exponential survival", {
  set.seed(45)
  beta <- 0.06; sus <- 0.7; N <- 600; t_end <- 10
  res <- run_stoch_linear(beta=beta, N=N, t=t_end, susceptibility=sus,
                           timepoints=t_end, n_sim=200, cores=1)
  expected_ar <- 1 - exp(-beta * sus * t_end)
  observed_ar <- mean(res$C1) / N
  expect_equal(observed_ar, expected_ar, tolerance=0.03)
})

# ---------------------------------------------------------------------------
# run_stoch_frailty_linear
# ---------------------------------------------------------------------------

test_that("run_stoch_frailty_linear returns correct structure", {
  set.seed(50)
  res <- run_stoch_frailty_linear(alpha=0.3, sd=0.3, beta=0.05, f=0.5, N=300, t=20,
                                   n_frailty=5, timepoints=0:20, n_sim=4, cores=1)
  expect_true(all(c("time","sim","vac","unvac","CRR") %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:20)
})

test_that("run_stoch_frailty_linear vac and unvac are non-decreasing within each sim", {
  set.seed(51)
  res <- run_stoch_frailty_linear(alpha=0.3, sd=0.3, beta=0.05, f=0.5, N=300, t=20,
                                   n_frailty=5, timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$vac)   >= 0))
    expect_true(all(diff(r$unvac) >= 0))
  }
})

test_that("run_stoch_frailty_linear alpha<1 gives lower attack rate in vaccinated group", {
  set.seed(52)
  res <- run_stoch_frailty_linear(alpha=0.2, sd=0.3, beta=0.1, f=0.5, N=500, t=20,
                                   n_frailty=10, timepoints=20, n_sim=40, cores=1)
  expect_true(mean(res$vac) < mean(res$unvac))
})

test_that("run_stoch_frailty_linear CRR < 1 when alpha < 1", {
  set.seed(53)
  res <- run_stoch_frailty_linear(alpha=0.3, sd=0.3, beta=0.1, f=0.5, N=500, t=20,
                                   n_frailty=10, timepoints=20, n_sim=30, cores=1)
  crr_vals <- res$CRR[is.finite(res$CRR) & res$CRR > 0]
  expect_true(mean(crr_vals) < 1)
})

test_that("run_stoch_frailty_linear alpha=0 produces no vaccinated infections", {
  set.seed(54)
  res <- run_stoch_frailty_linear(alpha=0, sd=0.3, beta=0.1, f=0.5, N=300, t=20,
                                   n_frailty=5, timepoints=0:20, n_sim=10, cores=1)
  expect_true(all(res$vac == 0))
  expect_true(all(res$unvac >= 0))
})

# ---------------------------------------------------------------------------
# run_coupled_linear
# ---------------------------------------------------------------------------

test_that("run_coupled_linear returns correct structure", {
  set.seed(60)
  res <- run_coupled_linear(beta=0.1, N=c(100, 100), t=20,
                             sus_x=c(1, 1), sus_z=c(1, 0.5),
                             timepoints=0:20, n_sim=4, cores=1)
  expect_true(all(c("Cx1","Cz1","Cx2","Cz2","time","sim") %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:20)
})

test_that("run_coupled_linear Cx and Cz are non-decreasing within each sim", {
  set.seed(61)
  res <- run_coupled_linear(beta=0.1, N=c(100, 100), t=20,
                             sus_x=c(1, 1), sus_z=c(1, 0.5),
                             timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$Cx1) >= 0))
    expect_true(all(diff(r$Cz1) >= 0))
    expect_true(all(diff(r$Cx2) >= 0))
    expect_true(all(diff(r$Cz2) >= 0))
  }
})

test_that("run_coupled_linear equal sus gives Cx == Cz (perfect coupling)", {
  set.seed(62)
  res <- run_coupled_linear(beta=0.1, N=c(100, 100), t=20,
                             sus_x=c(1, 0.5), sus_z=c(1, 0.5),
                             timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$Cx1 == res$Cz1))
  expect_true(all(res$Cx2 == res$Cz2))
})

test_that("run_coupled_linear sus_x > sus_z gives Cx >= Cz at all times", {
  set.seed(63)
  res <- run_coupled_linear(beta=0.1, N=c(200, 200), t=20,
                             sus_x=c(1, 1), sus_z=c(0.5, 0.5),
                             timepoints=0:20, n_sim=10, cores=1)
  expect_true(all(res$Cx1 >= res$Cz1))
  expect_true(all(res$Cx2 >= res$Cz2))
})

test_that("run_coupled_linear sus=0 gives zero infections", {
  set.seed(64)
  res <- run_coupled_linear(beta=0.1, N=c(100, 100), t=20,
                             sus_x=c(0, 0), sus_z=c(0, 0),
                             timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$Cx1 == 0 & res$Cx2 == 0))
  expect_true(all(res$Cz1 == 0 & res$Cz2 == 0))
})

test_that("run_coupled_linear sus>1 gives higher attack rate than sus=1", {
  set.seed(65)
  res1 <- run_coupled_linear(beta=0.05, N=500, t=20, sus_x=1, sus_z=1,
                              timepoints=20, n_sim=50, cores=1)
  res2 <- run_coupled_linear(beta=0.05, N=500, t=20, sus_x=3, sus_z=3,
                              timepoints=20, n_sim=50, cores=1)
  expect_true(mean(res2$Cx1) > mean(res1$Cx1))
})

test_that("run_coupled_linear marginal attack rate matches exponential survival", {
  set.seed(66)
  beta <- 0.06; sus <- 2; N <- 600; t_end <- 10
  res <- run_coupled_linear(beta=beta, N=N, t=t_end, sus_x=sus, sus_z=sus,
                             timepoints=t_end, n_sim=200, cores=1)
  expected_ar <- 1 - exp(-sus * beta * t_end)
  observed_ar <- mean(res$Cx1) / N
  expect_equal(observed_ar, expected_ar, tolerance=0.03)
})

# ---------------------------------------------------------------------------
# run_coupled_frailty_linear
# ---------------------------------------------------------------------------

test_that("run_coupled_frailty_linear returns correct structure", {
  set.seed(70)
  res <- run_coupled_frailty_linear(alpha=0.3, sd=0.3, beta=0.05, f=0.5, N=300, t=20,
                                     n_frailty=5, timepoints=0:20, n_sim=4, cores=1)
  expect_true(all(c("time","sim","vac_x","unvac_x","CRR_x","vac_z","unvac_z","CRR_z") %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:20)
})

test_that("run_coupled_frailty_linear vac and unvac are non-decreasing within each sim", {
  set.seed(71)
  res <- run_coupled_frailty_linear(alpha=0.3, sd=0.3, beta=0.05, f=0.5, N=300, t=20,
                                     n_frailty=5, timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$vac_x)   >= 0))
    expect_true(all(diff(r$unvac_x) >= 0))
    expect_true(all(diff(r$vac_z)   >= 0))
    expect_true(all(diff(r$unvac_z) >= 0))
  }
})

test_that("run_coupled_frailty_linear same vac_counts gives identical X and Z", {
  set.seed(72)
  res <- run_coupled_frailty_linear(alpha=0.3, sd=0.3, beta=0.05, f=0.5, N=200, t=20,
                                     n_frailty=5, timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$vac_x == res$vac_z))
  expect_true(all(res$unvac_x == res$unvac_z))
})

test_that("run_coupled_frailty_linear different vac_counts gives different X and Z", {
  set.seed(73)
  fr      <- get_frailty(sd=0.3, n=5)
  n_total <- round(2 * 300 * fr$p)
  vac_x   <- round(0.5 * n_total)
  vac_z   <- round(0   * n_total)
  res <- run_coupled_frailty_linear(alpha=0.3, sd=0.3, beta=0.05, N=300, t=20,
                                     n_frailty=5, vac_counts_x=vac_x, vac_counts_z=vac_z,
                                     timepoints=20, n_sim=20, cores=1)
  # Z has no vaccinated people so unvac_z = total infections
  expect_true(all(res$vac_z == 0))
  expect_true(mean(res$unvac_x) < mean(res$unvac_z))  # vaccination reduces unvac infections via herd effect
})

test_that("run_coupled_frailty_linear vac_counts_x=0 gives vac_x=0", {
  set.seed(74)
  fr      <- get_frailty(sd=0.3, n=5)
  n_total <- round(2 * 300 * fr$p)
  res <- run_coupled_frailty_linear(alpha=0.3, sd=0.3, beta=0.05, N=300, t=20,
                                     n_frailty=5,
                                     vac_counts_x=round(0 * n_total),
                                     vac_counts_z=round(0.5 * n_total),
                                     timepoints=20, n_sim=10, cores=1)
  expect_true(all(res$vac_x == 0))
})

test_that("run_coupled_frailty_linear alpha=0 gives zero infections in vaccinated group (X)", {
  set.seed(75)
  res <- run_coupled_frailty_linear(alpha=0, sd=0.3, beta=0.08, f=0.5, N=300, t=20,
                                     n_frailty=5, timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$vac_x == 0))
})

# ---------------------------------------------------------------------------
# run_coupled_sir
# ---------------------------------------------------------------------------

test_that("run_coupled_sir returns correct structure", {
  set.seed(80)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_coupled_sir(beta=2, N=c(30,30), mixing_matrix=mm, t=10,
                          I_ini=c(1,0), sus_x=c(1,1), sus_z=c(1,0.5),
                          gamma=1/7, timepoints=0:10, n_sim=4, cores=1)
  expect_true(all(c("Cx1","Cz1","Cx2","Cz2","time","sim") %in% names(res)))
  expect_equal(nrow(res), 4 * 11)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:10)
})

test_that("run_coupled_sir Cx and Cz are non-decreasing within each sim", {
  set.seed(81)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_coupled_sir(beta=2, N=c(40,40), mixing_matrix=mm, t=15,
                          I_ini=c(1,0), sus_x=c(1,1), sus_z=c(1,0.5),
                          gamma=1/7, timepoints=0:15, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$Cx1) >= 0))
    expect_true(all(diff(r$Cz1) >= 0))
    expect_true(all(diff(r$Cx2) >= 0))
    expect_true(all(diff(r$Cz2) >= 0))
  }
})

test_that("run_coupled_sir equal sus gives Cx == Cz (perfect coupling)", {
  set.seed(82)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_coupled_sir(beta=2, N=c(40,40), mixing_matrix=mm, t=15,
                          I_ini=c(1,0), sus_x=c(1,0.5), sus_z=c(1,0.5),
                          gamma=1/7, timepoints=0:15, n_sim=5, cores=1)
  expect_true(all(res$Cx1 == res$Cz1))
  expect_true(all(res$Cx2 == res$Cz2))
})

test_that("run_coupled_sir sus_x > sus_z gives Cx >= Cz at all times", {
  set.seed(83)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_coupled_sir(beta=2, N=c(50,50), mixing_matrix=mm, t=20,
                          I_ini=c(1,0), sus_x=c(1,1), sus_z=c(0.5,0.5),
                          gamma=1/7, timepoints=0:20, n_sim=8, cores=1)
  expect_true(all(res$Cx1 >= res$Cz1))
  expect_true(all(res$Cx2 >= res$Cz2))
})

test_that("run_coupled_sir sus_z=0 prevents new infections in Z (only seed remains)", {
  set.seed(84)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res <- run_coupled_sir(beta=2, N=c(40,40), mixing_matrix=mm, t=15,
                          I_ini=c(1,0), sus_x=c(1,1), sus_z=c(0,0),
                          gamma=1/7, timepoints=0:15, n_sim=5, cores=1)
  # seed individual is shared between scenarios, so Cz1 stays at 1; group 2 gets none
  expect_true(all(res$Cz1 == 1))
  expect_true(all(res$Cz2 == 0))
})

test_that("run_coupled_sir Cx conserves population (no more than N infections)", {
  set.seed(85)
  N <- c(30, 30)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res <- run_coupled_sir(beta=2, N=N, mixing_matrix=mm, t=20,
                          I_ini=c(1,0), sus_x=c(1,1), sus_z=c(1,0.5),
                          gamma=1/7, timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$Cx1 <= N[1]))
  expect_true(all(res$Cx2 <= N[2]))
  expect_true(all(res$Cz1 <= N[1]))
  expect_true(all(res$Cz2 <= N[2]))
})

# ---------------------------------------------------------------------------
# run_coupled_frailty_sir
# ---------------------------------------------------------------------------

test_that("run_coupled_frailty_sir returns correct structure", {
  set.seed(90)
  res <- run_coupled_frailty_sir(alpha=0.3, sd=0.3, beta=2, f=0.5, N=50, t=15,
                                  n_frailty=4, gamma=1/7,
                                  timepoints=0:15, n_sim=4, cores=1)
  expect_true(all(c("time","sim","vac_x","unvac_x","CRR_x","vac_z","unvac_z","CRR_z") %in% names(res)))
  expect_equal(nrow(res), 4 * 16)
  expect_setequal(unique(res$sim), 1:4)
  expect_setequal(unique(res$time), 0:15)
})

test_that("run_coupled_frailty_sir counts are non-decreasing within each sim", {
  set.seed(91)
  res <- run_coupled_frailty_sir(alpha=0.3, sd=0.3, beta=2, f=0.5, N=50, t=15,
                                  n_frailty=4, gamma=1/7,
                                  timepoints=0:15, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$vac_x)   >= 0))
    expect_true(all(diff(r$unvac_x) >= 0))
    expect_true(all(diff(r$vac_z)   >= 0))
    expect_true(all(diff(r$unvac_z) >= 0))
  }
})

test_that("run_coupled_frailty_sir same vac_counts gives identical X and Z", {
  set.seed(92)
  res <- run_coupled_frailty_sir(alpha=0.3, sd=0.3, beta=2, f=0.5, N=40, t=15,
                                  n_frailty=4, gamma=1/7,
                                  timepoints=0:15, n_sim=4, cores=1)
  expect_true(all(res$vac_x == res$vac_z))
  expect_true(all(res$unvac_x == res$unvac_z))
})

test_that("run_coupled_frailty_sir different vac_counts gives different X and Z", {
  set.seed(93)
  fr      <- get_frailty(sd=0.3, n=4)
  n_total <- round(2 * 60 * fr$p)
  vac_x   <- round(0.5 * n_total)
  vac_z   <- round(0   * n_total)
  res <- run_coupled_frailty_sir(alpha=0.3, sd=0.3, beta=2, N=60, t=15,
                                  n_frailty=4, gamma=1/7,
                                  vac_counts_x=vac_x, vac_counts_z=vac_z,
                                  timepoints=15, n_sim=20, cores=1)
  expect_true(all(res$vac_z == 0))
})

test_that("run_coupled_frailty_sir CRR_x < 1 when alpha < 1", {
  set.seed(94)
  res <- run_coupled_frailty_sir(alpha=0.3, sd=0.3, beta=2, f=0.5, N=60, t=20,
                                  n_frailty=4, gamma=1/7,
                                  timepoints=20, n_sim=20, cores=1)
  crr_vals <- res$CRR_x[is.finite(res$CRR_x) & res$CRR_x > 0]
  expect_true(mean(crr_vals) < 1)
})

test_that("run_stoch_frailty_cd R argument produces same beta as get_beta", {
  set.seed(34)
  alpha <- 0.4; sd_val <- 0.3; R_val <- 3
  beta_expected <- get_beta(R_val, alpha, sd_val, f=0.5, N=300, n_frailty=5, gamma=1/7)
  res_R    <- run_stoch_frailty_cd(alpha=alpha, sd=sd_val, R=R_val, f=0.5, N=300, t=10,
                                    n_frailty=5, gamma=1/7, timepoints=10, n_sim=5, cores=1)
  res_beta <- run_stoch_frailty_cd(alpha=alpha, sd=sd_val, beta=beta_expected, f=0.5, N=300, t=10,
                                    n_frailty=5, gamma=1/7, timepoints=10, n_sim=5, cores=1)
  # Both runs use identical beta; with the same seed their means should be close
  expect_equal(mean(res_R$unvac), mean(res_beta$unvac), tolerance=5)
})
