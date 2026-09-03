# Tests for the heterogeneous vaccine effect (alpha varying by individual):
# alpha_ve_bins(), run_stoch_ve_hetero_cd(), get_stoch_eate_ve_hetero().
#
# The headline question these were written for: does spreading alpha across
# individuals change the spread of the CIR ACROSS REALISATIONS, or only its
# level? It does both, downwards. Measured at N=400, t*=8, mean alpha=0.4,
# n_sim=6000: sd(CIR) 0.0899 -> 0.0749 (-17%) as kappa goes 0 -> 0.8.
#
# Why it falls, since "more heterogeneity" sounds like it should mean more
# variability. Under homogeneous mixing both arms see one common cumulative
# FOI Lambda, and
#     R(Lambda) = [1 - sum_k p_k exp(-alpha_k Lambda)] / [1 - exp(-Lambda)]
# Two mechanisms, both pushing the same way:
#
#  1. R(Lambda) FLATTENS. R normally rises with Lambda (both arms saturate
#     towards 1 -- the familiar "VE looks worse in bigger epidemics"
#     artefact). Spreading alpha cuts dR/dLambda from 0.122 to 0.044, a 64%
#     flattening, because the epidemic selects: the poorly protected
#     vaccinated are infected first, so the surviving vaccinated pool is
#     progressively enriched for the well protected, buffering that arm
#     against epidemic size. Homogeneity has nothing to select on.
#
#  2. Heterogeneous Bernoullis are UNDER-dispersed. At fixed Lambda,
#     Var(C2) = sum_k n_k p_k (1-p_k) = n[pbar(1-pbar) - Var_k(p_k)],
#     and spreading alpha spreads p_k, so that term is subtracted.
#
# Decomposing Var(CIR) by epidemic-size decile says which dominates: the
# between-realisation (epidemic size) component collapses 86%, but it is a
# small share of the total to begin with, so it is only ~23% of the variance
# drop; the within component supplies ~77%. sd(Lambda) itself moves just -9%.
#
# Caveat for anyone tempted to use sd(CIR) as a heterogeneity detector: the
# mean CIR falls too (0.488 -> 0.441), so about half the raw-sd effect is
# scale, not extra regularity -- the CV moves only 8%. The slope of CIR
# against attack rate is the much sharper signature (64% vs 17%), and being a
# shape effect, no rescaling mimics it.

skip_if_not_installed("data.table")

library(data.table)

if (!exists("run_stoch_ve_hetero_cd")) {
  source("../../utils.R", chdir = TRUE)
  source("../../stoch_model.R", chdir = TRUE)
}

# Shared scenario. Every simulator call passes an explicit dust seed, so these
# tests are deterministic rather than merely probable-to-pass.
VH <- list(N = 400, t_star = 8, gamma = 1, beta = 2.0, alpha = 0.4,
           n_alpha = 10L, I_ini_total = 10, dt = 0.01, seed = 11L)

vh_run <- function(kappa, n_sim, seed = VH$seed, n_alpha = VH$n_alpha) {
  run_stoch_ve_hetero_cd(
    beta = VH$beta, susceptibility = c(1, VH$alpha),
    alpha_kappa = kappa, n_alpha = n_alpha,
    N = VH$N / 2, t = VH$t_star, gamma = VH$gamma,
    I_ini_total = VH$I_ini_total, timepoints = seq(1, VH$t_star, 1),
    n_sim = n_sim, cores = 1, dt = VH$dt, f = 0.5, seed = seed)
}

# CIR = (C2 / N_vac) / (C1 / N_unvac), per realisation. It is a ratio of
# counts, so realisations where the epidemic never took off contribute 0/0 or
# a denominator near zero, and the sd is then set by a handful of exploding
# ratios rather than by the biology. Restrict to takeoff, and report how much
# was dropped so a silent change in that fraction cannot masquerade as a
# change in the sd.
vh_cir <- function(dt_out, fizzle_frac = 0.05) {
  f    <- dt_out[time == VH$t_star]
  keep <- (f$C1 + f$C2) > fizzle_frac * VH$N & f$C1 > 0
  list(cir = (f$C2[keep] / (VH$N / 2)) / (f$C1[keep] / (VH$N / 2)),
       frac_kept = mean(keep))
}

# ---------------------------------------------------------------------------
# alpha_ve_bins: the bin construction
# ---------------------------------------------------------------------------

test_that("bin mean equals alpha exactly, across alpha x kappa x n_alpha", {
  # The mean is the quantity the optimiser identifies. If it drifts, the fit
  # does not fail -- it converges on the wrong alpha, which is far worse.
  # Accumulate over the grid and assert once, rather than emitting a few
  # thousand expectations.
  worst_mean <- 0; bad_len <- 0L; bad_p <- 0; out_of_range <- 0L
  for (K in c(1L, 5L, 10L, 20L)) {
    for (kappa in c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95)) {
      for (a in seq(0.02, 0.98, 0.04)) {
        b <- alpha_ve_bins(a, kappa, K)
        if (length(b$x) != K) bad_len <- bad_len + 1L
        bad_p        <- max(bad_p, abs(sum(b$p) - 1))
        out_of_range <- out_of_range + sum(b$x < 0 | b$x > 1)
        worst_mean   <- max(worst_mean, abs(sum(b$p * b$x) - a))
      }
    }
  }
  expect_equal(bad_len, 0L)
  expect_equal(out_of_range, 0L)
  expect_lt(bad_p, 1e-12)
  expect_lt(worst_mean, 1e-10)
})

test_that("kappa = 0 is a point mass, and spread grows with kappa", {
  b0 <- alpha_ve_bins(0.4, 0, 10L)
  expect_true(all(b0$x == 0.4))

  sds <- sapply(c(0, 0.1, 0.3, 0.5, 0.7),
                function(k) { b <- alpha_ve_bins(0.4, k, 20L)
                              sqrt(sum(b$p * (b$x - 0.4)^2)) })
  expect_true(all(diff(sds) > 0))
  # Target is sd = kappa * sqrt(a(1-a)); binning loses a little, but not much.
  target <- c(0, 0.1, 0.3, 0.5, 0.7) * sqrt(0.4 * 0.6)
  expect_true(all(sds <= target + 1e-9))
  expect_true(all(sds[-1] > 0.95 * target[-1]))
})

test_that("alpha_ve_bins rejects an out-of-range kappa", {
  expect_error(alpha_ve_bins(0.4, 1.0, 10L))
  expect_error(alpha_ve_bins(0.4, -0.1, 10L))
})

# ---------------------------------------------------------------------------
# kappa = 0 must be the homogeneous SIR
# ---------------------------------------------------------------------------

test_that("kappa = 0 with one bin reproduces the homogeneous SIR exactly", {
  # With n_alpha = 1 the group layout is literally the 2-group SIR, so this is
  # an exact-equality check, not a statistical one.
  ref <- run_stoch_cd_dust(matrix(1, 2, 2), beta = VH$beta,
                           N = c(VH$N / 2, VH$N / 2), t = VH$t_star,
                           I_ini = c(5, 5), susceptibility = c(1, VH$alpha),
                           gamma = VH$gamma, dt = VH$dt,
                           timepoints = seq(1, VH$t_star, 1),
                           n_sim = 300, cores = 1, seed = VH$seed)
  setDT(ref)
  got <- vh_run(kappa = 0, n_sim = 300, n_alpha = 1L)
  expect_equal(got[time == VH$t_star]$C1, ref[time == VH$t_star]$C1)
  expect_equal(got[time == VH$t_star]$C2, ref[time == VH$t_star]$C2)
})

# ---------------------------------------------------------------------------
# CIR variance across realisations -- the question this file exists for
# ---------------------------------------------------------------------------

test_that("spreading alpha lowers both the level and the spread of the CIR", {
  n_sim <- 2500
  ks    <- c(0, 0.4, 0.8)
  out   <- lapply(ks, function(k) vh_cir(vh_run(k, n_sim)))

  # Takeoff fraction must be stable across kappa, or a change in sd could just
  # be a change in which realisations survived the filter.
  fk <- vapply(out, function(o) o$frac_kept, numeric(1))
  expect_true(all(abs(fk - fk[1]) < 0.02))
  expect_true(all(vapply(out, function(o) length(o$cir), integer(1)) > 0.9 * n_sim))

  mean_cir <- vapply(out, function(o) mean(o$cir), numeric(1))
  sd_cir   <- vapply(out, function(o) sd(o$cir),   numeric(1))

  # Level: heterogeneity makes the vaccine look BETTER than its mean alpha.
  # exp(-a * Lambda) is convex in a, so E[exp(-a_i L)] > exp(-mean(a) L) --
  # more vaccinated survivors, fewer vaccinated cases, lower CIR. This is
  # Jensen, not noise, so it is a strict inequality.
  expect_true(all(diff(mean_cir) < 0))

  # Spread: sd(CIR) across realisations FALLS as alpha is spread out.
  # Measured at N=400, t*=8, mean alpha=0.4, n_sim=6000: sd 0.0899 -> 0.0749
  # from kappa 0 to 0.8, a 17% drop, with non-overlapping bootstrap CIs.
  # se(sd) ~ sd/sqrt(2n) is under 2% at this n_sim, so 17% is ~9 se -- the
  # 8% threshold below is a wide margin, not a tuned one.
  expect_lt(sd_cir[3], sd_cir[1])
  expect_lt(sd_cir[3] / sd_cir[1], 0.92)

  # Part of that drop is just the lower mean. The CV falls too, so it is not
  # ONLY a scale effect -- but the CV effect is roughly half the size, which
  # is why the raw-sd claim on its own would overstate it.
  cv <- sd_cir / mean_cir
  expect_lt(cv[3], cv[1])
})

test_that("the CIR sd effect survives a change of seed", {
  # Guards against the headline result being one lucky dust stream.
  for (sd_seed in c(101L, 202L)) {
    s0 <- sd(vh_cir(vh_run(0,   1500, seed = sd_seed))$cir)
    s8 <- sd(vh_cir(vh_run(0.8, 1500, seed = sd_seed))$cir)
    expect_lt(s8, s0)
  }
})

# ---------------------------------------------------------------------------
# EATE
# ---------------------------------------------------------------------------

test_that("EATE at kappa = 0 matches the homogeneous SIR EATE", {
  tp  <- seq(1, VH$t_star, 1)
  ref <- get_stoch_eate_sir(beta = VH$beta, susceptibility = c(1, VH$alpha),
                            f = 0.5, N = VH$N, t = VH$t_star, gamma = VH$gamma,
                            I_ini = c(5, 5), n_vac = 6, n_rep = 40,
                            dt = VH$dt, timepoints = tp, mc.cores = 1)
  got <- get_stoch_eate_ve_hetero(alpha = VH$alpha, alpha_kappa = 0,
                                  n_alpha = VH$n_alpha, beta = VH$beta,
                                  f = 0.5, N = VH$N / 2, t = VH$t_star,
                                  gamma = VH$gamma, I_ini_total = VH$I_ini_total,
                                  n_vac = 6, n_rep = 40, dt = VH$dt,
                                  timepoints = tp, mc.cores = 1)
  setDT(ref); setDT(got)
  ve <- function(d) 1 - mean(d[method == "full_stoch" & t == VH$t_star]$eate)
  expect_equal(ve(got), ve(ref), tolerance = 0.05)
})

test_that("EATE VE rises with kappa at fixed mean alpha", {
  # The counterfactual is per bin: a flipped individual gets THEIR alpha_k.
  # Collapsing to a mean alpha first would flatten this out entirely, so a
  # rising VE is what shows the per-bin counterfactual is actually wired up.
  tp <- seq(1, VH$t_star, 1)
  ve <- vapply(c(0, 0.8), function(k) {
    g <- get_stoch_eate_ve_hetero(alpha = VH$alpha, alpha_kappa = k,
                                  n_alpha = VH$n_alpha, beta = VH$beta,
                                  f = 0.5, N = VH$N / 2, t = VH$t_star,
                                  gamma = VH$gamma, I_ini_total = VH$I_ini_total,
                                  n_vac = 8, n_rep = 50, dt = VH$dt,
                                  timepoints = tp, mc.cores = 1)
    setDT(g)
    1 - mean(g[method == "full_stoch" & t == VH$t_star]$eate)
  }, numeric(1))
  expect_gt(ve[2], ve[1] + 0.02)
})
