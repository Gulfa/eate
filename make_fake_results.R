# Synthetic fit-array results, for exercising analyse_fit_array.R locally
# ---------------------------------------------------------------------------
# The analysis script is hard to test without a cluster run: it wants a
# directory of results_*.RDS whose jobs carry a dozen interlocking fields, and a
# layout change can break it in ways that only surface at the very end. This
# builds a small set of jobs in exactly that shape -- linear, two SIR I_ini
# variants, three network jobs across two pl_alphas, multisite and parity --
# with planted VE trends over t and over coverage, so every panel has something
# to draw.
#
#   Rscript make_fake_results.R                       # writes output/_fake_res
#   Rscript analyse_fit_array.R output/_fake_res      # then run the analysis
#
# The numbers are arbitrary; only the SHAPE is meaningful. Never read a result
# off this.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages(library(data.table)))
set.seed(21)
tps <- c(0.5, 1:8); t_star <- 8; K <- 8; S <- 3

mkrows <- function(base, covv = NULL) {
  d <- rbindlist(lapply(1:K, function(k) {
    rbindlist(lapply(tps, function(tt) {
      data.table(
        method = rep(c("full_stoch", "CRR"), each = S),
        t = tt, param_sample = k, sim = rep(1:S, 2),
        eate = c(1 - (base - 0.012 * tt + rnorm(S, 0, .02)),
                 1 - (base - 0.02 + rnorm(S, 0, .02))),
        ave  = c(0.18 - 0.01 * tt + rnorm(S, 0, .01), rep(0.1, S)),
        num = 100, denom = 200,
        beta_k = rnorm(1, 2.1, 0.1), alpha_k = rnorm(1, 1 - base, 0.02),
        eate_sd_rep = c(rep(0.15, S), rep(NA_real_, S)),
        ave_sd_rep  = c(rep(0.06, S), rep(NA_real_, S)))
    }))
  }))
  if (!is.null(covv)) d[, coverage := covv]
  d[]
}

mkjob <- function(nm, mt, pa = NULL, base = 0.5, net_seed = NULL, alloc = 1L,
                  chisq = 1.2) {
  list(name = nm, experiment_id = "expT", model_type = mt, fit_method = "kernel",
       pl_alpha = pa, network_seed = net_seed, allocation_seed = alloc,
       network_engine = if (identical(mt, "network")) "events" else NA_character_,
       ve_cf_method = "resim",
       ve_n_flip_used = if (identical(mt, "network")) 100L else NA_integer_,
       parity_alpha_alt = NA_real_, parity_mod = NA_integer_,
       parity_alpha_up = NA_real_, parity_alpha_down = NA_real_,
       split_frac = NA_real_, split_alpha_prod = NA_real_,
       alpha_kappa = NA_real_, n_alpha = NA_integer_,
       fit = list(beta = 2.1, alpha = 1 - base, loss = 0.5, convergence = 0L),
       fit_bounds = list(beta = c(0.05, 5), alpha = c(0.01, 2)),
       posterior_cov = list(cov = diag(2), J = NULL, Sigma = NULL,
                            sd = c(beta = 0.1, alpha = 0.02)),
       grid_post = NULL, design_coverage = 0.5,
       loss_floor = 0.4, loss_chisq = chisq, resid_C1 = 2, resid_C2 = -1,
       ve = { d <- mkrows(base); d[, model := mt]; d[] },
       ve_uncertainty = mkrows(base),
       ve_by_coverage = rbindlist(lapply(c(0.25, 0.75), function(f)
         mkrows(base + 0.05 * (f - 0.5), f))),
       coverage_effect = data.table(averted_per1k = rnorm(K, 20, 3),
                                    cov_from = 0.5, cov_to = 0.6))
}

jobs <- list(
  mkjob("expT__linear",    "linear",        base = 0.44),
  mkjob("expT__sir_i10",   "sir_i10",       base = 0.41),
  mkjob("expT__sir_i20",   "sir_i20",       base = 0.42),
  mkjob("expT__net14_n1",  "network", pa = 1.4, base = 0.52, net_seed = 1L),
  mkjob("expT__net14_n2",  "network", pa = 1.4, base = 0.51, net_seed = 2L),
  mkjob("expT__net3_n1",   "network", pa = 3,   base = 0.47, net_seed = 1L),
  mkjob("expT__net5_bad",  "network", pa = 5,   base = 0.60, net_seed = 1L,
        chisq = 14.2),
  mkjob("expT__ms",        "sir_multisite", base = 0.45),
  mkjob("expT__ms_bad",    "sir_multisite", base = 0.66, alloc = 2L, chisq = 9.1),
  mkjob("expT__par",       "sir_parity_m2", base = 0.30))

dir.create("output/_fake_res", recursive = TRUE, showWarnings = FALSE)
saveRDS(jobs, "output/_fake_res/results_01.RDS")
cat("wrote", length(jobs), "synthetic jobs\n")
