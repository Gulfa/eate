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
# run_stoch_cd
# ---------------------------------------------------------------------------

test_that("run_stoch_cd returns correct structure for n groups", {
  set.seed(20)
  mm <- matrix(1/4, nrow=4, ncol=4)
  res <- run_stoch_cd(mm, beta=1, N=c(100,100,100,100), t=20, I_ini=c(1,0,0,0),
                      timepoints=0:20, n_sim=4, cores=1)
  expected_cols <- c(paste0(c("S","I","R","C"), rep(1:4, each=1)) |>
                       (\(x) c(paste0("S",1:4), paste0("I",1:4), paste0("R",1:4), paste0("C",1:4)))(),
                     "time", "sim")
  expect_true(all(expected_cols %in% names(res)))
  expect_equal(nrow(res), 4 * 21)
  expect_setequal(unique(res$sim), 1:4)
})

test_that("run_stoch_cd conserves population in every group", {
  set.seed(21)
  mm <- matrix(c(0.8,0.2, 0.2,0.8), nrow=2)
  N  <- c(150, 100)
  res <- run_stoch_cd(mm, beta=1.5, N=N, t=30, I_ini=c(2,0),
                      timepoints=0:30, n_sim=10, cores=1)
  expect_true(all(res$S1 + res$I1 + res$R1 == N[1]))
  expect_true(all(res$S2 + res$I2 + res$R2 == N[2]))
})

test_that("run_stoch_cd cumulative cases are non-decreasing within each sim", {
  set.seed(22)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd(mm, beta=1.5, N=c(200,200), t=20, I_ini=c(2,0),
                      timepoints=0:20, n_sim=5, cores=1)
  for (s in unique(res$sim)) {
    r <- res[res$sim == s, ]
    expect_true(all(diff(r$C1) >= 0))
    expect_true(all(diff(r$C2) >= 0))
  }
})

test_that("run_stoch_cd beta=0 produces no infections", {
  set.seed(23)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd(mm, beta=0, N=c(100,100), t=20, I_ini=c(0,0),
                      timepoints=0:20, n_sim=5, cores=1)
  expect_true(all(res$C1 == 0))
  expect_true(all(res$C2 == 0))
})

test_that("run_stoch_cd susceptibility=0 for one group blocks its infections", {
  set.seed(24)
  mm <- matrix(1/2, nrow=2, ncol=2)
  res <- run_stoch_cd(mm, beta=2, N=c(200,200), t=30, I_ini=c(3,0),
                      susceptibility=c(1, 0), timepoints=30, n_sim=20, cores=1)
  expect_true(all(res$C2 == 0))
  expect_true(mean(res$C1) > 10)
})

test_that("run_stoch_cd with homogeneous mixing and equal susceptibility matches run_sir attack rate", {
  set.seed(25)
  mm  <- matrix(1/2, nrow=2, ncol=2)
  res_cd  <- run_stoch_cd(mm, beta=1.5, N=c(200,200), t=40, I_ini=c(2,2),
                           susceptibility=c(1,1), gamma=1/7,
                           timepoints=40, n_sim=80, cores=1)
  res_sir <- run_sir(N_cont=200, N_vac=200, beta=1.5, alpha=1, gamma=1/7,
                     I0_cont=2, I0_vac=2, t=40, timepoints=40, n_sim=80, cores=1)
  ar_cd  <- (mean(res_cd$C1)  + mean(res_cd$C2))  / 400
  ar_sir <- (mean(res_sir$C1) + mean(res_sir$C2)) / 400
  expect_equal(ar_cd, ar_sir, tolerance=0.05)
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
