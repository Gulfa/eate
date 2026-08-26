# Slurm array runner that combines optim_histograms.R (fit + posterior
# cov) with network_alpha_sweep.R (vary network_seed and allocation_seed).
#
# For each job (one config), do:
#   1. Build the per-model stochastic simulator (linear / sir /
#      sus-frailty / trans-frailty / network), with c_ij and vac
#      allocation materialised from the (network_seed, allocation_seed)
#      embedded in the config.
#   2. Grid-search + Nelder-Mead (with restarts) to fit (beta, alpha)
#      so that mean(C1), mean(C2) at t* match (data_C1, data_C2).
#   3. Estimate posterior covariance via J^-1 Sigma J^-T with central
#      finite differences and CRN (shared dust seed across all 5
#      simulator calls).
#   4. Compute VE(t) at the fitted (beta, alpha) using the matching
#      EATE function (deterministic for linear/sir, stochastic frozen-
#      field for frailty/network).
#
# Submit with e.g.
#   sbatch --array=1-20 --cpus-per-task=4 run_fit_array.sh
# where run_fit_array.sh wraps `Rscript run_fit_array.R`.

library(dplyr)
library(data.table)
library(glue)

source("det_model.R")
source("stoch_model.R")
source("utils.R")

RhpcBLASctl::omp_set_num_threads(1)
setDTthreads(1)

# ---------------------------------------------------------------------------
# Experiments and shared knobs
# ---------------------------------------------------------------------------
#
# Each "experiment" is a scenario with its own population, targets, and
# initial-infection seeds. Every experiment gets its own set of configs
# spanning all model types (linear / sir / frailty / network) — nothing
# is pooled across experiments in the plots.

experiments <- list(
  list(id = "expA_N200_t8",
       N_cont = 100, N_vac = 100,
       data_C1 = 50, data_C2 = 25,
       t_star = 8,
       I_ini_2g = c(10, 10), init_I_nw = 2),
  list(id = "expB_N400_t12",
       N_cont = 200, N_vac = 200,
       data_C1 = 100, data_C2 = 50,
       t_star = 12,
       I_ini_2g = c(20, 20), init_I_nw = 4)
)

gamma  <- 1
dt     <- 0.01

# Optim
n_sim_opt   <- 1000
optim_maxit <- 250
n_restarts  <- 3
restart_loss_threshold <- 5
grid_n      <- 6
# Common Random Numbers for the fit loss: every loss evaluation uses the
# SAME fixed set of seeds, so the loss is a deterministic, smooth function
# of (beta, alpha) rather than being re-randomised at each optim probe.
# This is what makes Nelder-Mead converge reproducibly (mirrors the CRN
# already used in estimate_posterior_cov). Averaging over n_seed_opt seeds
# smooths the surface and reduces the seed-specific bias in (beta, alpha).
opt_seed    <- 20240517L
n_seed_opt  <- 3
log_beta_lo  <- log(0.01); log_beta_hi  <- log(5)
log_alpha_lo <- log(0.01); log_alpha_hi <- log(2)

# Fit method. "mean": least-squares on the simulated mean (correct when the
# final-size distribution is unimodal; cheapest). "kernel": negative-log kernel
# synthetic likelihood, which matches the data to the MODE and is unbiased when
# the distribution is bimodal (low initial infecteds, e.g. the network's
# init_I_nw = 2). The kernel posterior comes from the kernel-Hessian with an
# h->0 extrapolation and needs a larger n_sim (the MC-likelihood Hessian is
# noisy). Assumes data_C1/C2 are a single observed realisation, not E[C].
fit_method       <- "mean"
kernel_h         <- 3            # kernel bandwidth (count units) for the fit loss
kernel_hess_hs   <- c(2, 3, 4)   # bandwidths for the h->0 posterior extrapolation
kernel_hess_nsim <- 8000L        # sims for the (noisy) kernel Hessian
kernel_hess_hrel <- 0.03         # FD step (relative) for the 2nd-derivative Hessian

# Posterior cov
post_cov_n_sim <- 1000
post_cov_seed  <- 1234L
post_cov_h     <- 0.01

# VE (stochastic EATE)
ve_n_vac <- 5
ve_n_rep <- 200

# Posterior uncertainty propagation: K parameter samples drawn from
# MVN(fit, posterior_cov), each fed through compute_ve. Use a smaller
# n_rep on this pass since K * n_vac inner allocations already average
# out a lot of per-rep noise.
K_post_samples   <- 30
ve_n_rep_uncert  <- 50

# Parallelism — only two knobs; the code decides how to spend them
# per phase (see run_one_job below).
#
#   cores_per_node : total CPU budget per slurm task (--cpus-per-task).
#   nodes          : number of slurm array tasks (--array=1-nodes).
#
# Per-phase strategy:
#   Fit + posterior cov + point VE : all dust calls are sequential in R,
#     so `inner_cores = cores_per_node` gives each dust call the full
#     thread budget.
#   VE with uncertainty            : dominant sequential K loop; run K
#     posterior samples in parallel via mclapply with `K_cores =
#     cores_per_node`, and inside each K worker use `inner_cores = 1`
#     (each sim call single-threaded, no over-subscription).
#
# Configs within a slurm task run SEQUENTIALLY (no outer mclapply):
# parallelism happens inside each config; cross-config parallelism is
# via slurm array tasks (nodes).
cores_per_node <- 32
nodes          <- 20

out_dir <- "output/fit_array_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Build the config grid:
#   - one config each for linear, sir, sir_sus_frailty, sir_trans_frailty
#   - n_networks x n_allocations configs for network
# ---------------------------------------------------------------------------

n_networks              <- 10
n_allocations           <- 10   # outer allocations per network_seed
n_allocations_frailty   <- 10   # outer allocations per frailty config
n_allocations_multisite <- 10   # outer allocations per multi-site config
pl_alphas               <- c(2, 3, 5)   # Pareto exponents to sweep (kept separate)
mean_k                  <- 6

# Multi-site RCT knobs. n_sites locations; site_icc in [0,1] is the intra-site
# correlation of vaccine status. 0 = simple (individual) randomisation:
# per-site vaccinated counts vary by chance around f with the global total
# held exact (the realistic phase-3 design); 1 = cluster randomised (whole
# sites single-arm). Edit to vary.
multisite_n_sites <- 4
multisite_icc     <- 0

base_common <- list(
  gamma = gamma, dt = dt,
  n_sim_opt = n_sim_opt, optim_maxit = optim_maxit,
  n_restarts = n_restarts, restart_loss_threshold = restart_loss_threshold,
  grid_n = grid_n, opt_seed = opt_seed, n_seed_opt = n_seed_opt,
  # inner_cores is set by run_one_job per phase from cores_per_node.
  inner_cores = 1L,
  post_cov_n_sim = post_cov_n_sim, post_cov_seed = post_cov_seed,
  post_cov_h = post_cov_h,
  fit_method = fit_method, kernel_h = kernel_h,
  kernel_hess_hs = kernel_hess_hs, kernel_hess_nsim = kernel_hess_nsim,
  kernel_hess_hrel = kernel_hess_hrel,
  ve_n_vac = ve_n_vac, ve_n_rep = ve_n_rep,
  K_post_samples = K_post_samples, ve_n_rep_uncert = ve_n_rep_uncert
)

build_configs_for_experiment <- function(exp) {
  base <- modifyList(base_common,
                     list(experiment_id = exp$id,
                          N_cont = exp$N_cont, N_vac = exp$N_vac,
                          data_C1 = exp$data_C1, data_C2 = exp$data_C2,
                          t_star = exp$t_star,
                          I_ini_2g = exp$I_ini_2g, init_I_nw = exp$init_I_nw))
  cs <- list()

  # Linear and homogeneous SIR: allocation has no effect (exchangeable
  # individuals), so a single config each + ve_n_vac = 1 inside the EATE
  # function.
  cs[[length(cs)+1]] <- modifyList(base, list(
    name = glue("{exp$id}__linear"),
    model_type = "linear", ve_n_vac = 1))
  cs[[length(cs)+1]] <- modifyList(base, list(
    name = glue("{exp$id}__sir"),
    model_type = "sir", ve_n_vac = 1))
  # Multi-site RCT: n_sites locations, per-site vaccine fraction dispersion
  # set by site_icc (0 = individually randomised / max within-site
  # cancellation, 1 = cluster randomised / fully separated). The allocation
  # (which sites cluster how) varies by seed, so loop it like the frailty
  # models. ve_n_vac inherits base (draws fresh site allocations for VE).
  for (alloc_seed in seq_len(n_allocations_multisite)) {
    cs[[length(cs)+1]] <- modifyList(base, list(
      name            = glue("{exp$id}__sir_multisite_a{alloc_seed}"),
      model_type      = "sir_multisite",
      n_sites = multisite_n_sites, site_icc = multisite_icc,
      allocation_seed = alloc_seed))
  }

  # Frailty models: allocation matters (which bins get vaccinated).
  # Per-bin susceptibility is exp(frailty_amp * x), x ~ Beta on [0,1]
  # (mean-1 normalised). The Beta `sd` (capped < 0.5) sets how bimodal x is;
  # `frailty_amp` sets how strongly that maps into susceptibility, so it is
  # the unbounded lever for total frailty heterogeneity.
  for (alloc_seed in seq_len(n_allocations_frailty)) {
    cs[[length(cs)+1]] <- modifyList(base, list(
      name            = glue("{exp$id}__sir_sus_frailty_a{alloc_seed}"),
      model_type      = "sir_sus_frailty",
      sd = 0.3, sd_trans = 0, n_frailty = 10, frailty_amp = 2.5,
      allocation_seed = alloc_seed))
    cs[[length(cs)+1]] <- modifyList(base, list(
      name            = glue("{exp$id}__sir_trans_frailty_a{alloc_seed}"),
      model_type      = "sir_trans_frailty",
      sd = 0, sd_trans = 0.3, n_frailty = 10, frailty_amp = 2.5,
      allocation_seed = alloc_seed))
  }

  # Network: sweep pl_alpha as an outer axis, kept separate in all plots.
  for (pa in pl_alphas) {
    for (network_seed in seq_len(n_networks)) {
      for (alloc_seed in seq_len(n_allocations)) {
        cs[[length(cs)+1]] <- modifyList(base, list(
          name            = glue("{exp$id}__network_pa{pa}_n{network_seed}_a{alloc_seed}"),
          model_type      = "network",
          pl_alpha        = pa, mean_k = mean_k,
          network_seed    = network_seed,
          allocation_seed = alloc_seed))
      }
    }
  }
  cs
}

configs <- unlist(lapply(experiments, build_configs_for_experiment),
                  recursive = FALSE)

message(glue("Built {length(configs)} configs."))

# ---------------------------------------------------------------------------
# Per-model simulator builder. Returns function(beta, alpha, n_sim, seed)
# producing a data.table at time == t_star with columns (sim, C1, C2).
# ---------------------------------------------------------------------------

at_tstar <- function(out, t_star) {
  setDT(out)
  out[time == t_star, .(sim, C1, C2)]
}

build_simulator <- function(cfg) {
  N_total <- cfg$N_cont + cfg$N_vac
  switch(cfg$model_type,
    linear = function(beta, alpha, n_sim, seed = NULL) {
      at_tstar(run_stoch_linear_dust(
        beta = beta, N = c(cfg$N_cont, cfg$N_vac),
        susceptibility = c(1, alpha),
        t = cfg$t_star, dt = cfg$dt,
        timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores, seed = seed), cfg$t_star)
    },
    sir = function(beta, alpha, n_sim, seed = NULL) {
      at_tstar(run_stoch_cd_dust(
        matrix(rep(1, 4), nrow = 2),
        beta = beta, N = c(cfg$N_cont, cfg$N_vac),
        t = cfg$t_star, I_ini = cfg$I_ini_2g,
        susceptibility = c(1, alpha),
        gamma = cfg$gamma, dt = cfg$dt,
        timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores, seed = seed), cfg$t_star)
    },
    sir_multisite = function(beta, alpha, n_sim, seed = NULL) {
      # L sites, block-diagonal; per-site vaccine allocation from
      # cfg$.vac_sites (materialised from site_icc). Pool cases across sites:
      # C1 = sum unvaccinated cases, C2 = sum vaccinated cases. site_icc
      # controls within-site cancellation (0 = balanced, 1 = clustered).
      L          <- cfg$n_sites
      N_site_vec <- multisite_site_sizes(cfg$N_cont + cfg$N_vac, L)
      vac_l      <- cfg$.vac_sites
      N_grp      <- as.integer(rbind(N_site_vec - vac_l, vac_l))
      sus_grp    <- as.numeric(rbind(rep(1, L), rep(alpha, L)))
      mm         <- multisite_block_matrix(L, L)
      I_grp      <- .spread_seeds(sum(cfg$I_ini_2g), N_grp)
      out <- run_stoch_cd_dust(
        mm, beta = beta, N = N_grp, t = cfg$t_star, I_ini = I_grp,
        susceptibility = sus_grp, gamma = cfg$gamma, dt = cfg$dt,
        timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores, seed = seed)
      setDT(out)
      fin <- out[time == cfg$t_star]
      uc <- paste0("C", seq(1, 2L * L, 2)); vc <- paste0("C", seq(2, 2L * L, 2))
      data.table(sim = fin$sim,
                 C1 = rowSums(as.matrix(fin[, ..uc])),
                 C2 = rowSums(as.matrix(fin[, ..vc])))
    },
    sir_sus_frailty = function(beta, alpha, n_sim, seed = NULL) {
      at_tstar(run_stoch_frailty_cd(
        sd = cfg$sd, sd_trans = cfg$sd_trans,
        susceptibility = c(1, alpha), beta = beta,
        N = N_total / 2, t = cfg$t_star,
        n_frailty = cfg$n_frailty, gamma = cfg$gamma,
        I_ini_total = sum(cfg$I_ini_2g),
        vac_counts = cfg$.vac_counts,
        timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores,
        method = "dust", dt = cfg$dt, f = 0.5,
        frailty_amp = cfg$frailty_amp %||% 2.5, seed = seed), cfg$t_star)
    },
    sir_trans_frailty = function(beta, alpha, n_sim, seed = NULL) {
      at_tstar(run_stoch_frailty_cd(
        sd = cfg$sd, sd_trans = cfg$sd_trans,
        susceptibility = c(1, alpha), beta = beta,
        N = N_total / 2, t = cfg$t_star,
        n_frailty = cfg$n_frailty, gamma = cfg$gamma,
        I_ini_total = sum(cfg$I_ini_2g),
        vac_counts = cfg$.vac_counts,
        timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores,
        method = "dust", dt = cfg$dt, f = 0.5,
        frailty_amp = cfg$frailty_amp %||% 2.5, seed = seed), cfg$t_star)
    },
    network = function(beta, alpha, n_sim, seed = NULL) {
      at_tstar(run_stoch_network(
        beta = beta, N = N_total,
        susceptibility = c(1, alpha),
        t = cfg$t_star, c_ij = cfg$.c_ij, vac = cfg$.vac,
        k_mean = cfg$mean_k, gamma = cfg$gamma,
        dt = cfg$dt, timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores,
        I_ini = cfg$init_I_nw, seed = seed), cfg$t_star)
    },
    stop("Unknown model_type: ", cfg$model_type))
}

# Materialise per-config side state (network c_ij + vac allocation) so
# everything else can be configs of scalars.
materialise_cfg <- function(cfg) {
  if (cfg$model_type == "network") {
    set.seed(cfg$network_seed)
    cfg$.c_ij <- get_conact_matrix_pl(cfg$N_cont + cfg$N_vac,
                                       alpha = cfg$pl_alpha,
                                       mean_k = cfg$mean_k)
    set.seed(NULL)
    set.seed(cfg$allocation_seed)
    cfg$.vac <- sample(seq_len(cfg$N_cont + cfg$N_vac), cfg$N_vac)
    set.seed(NULL)
  } else if (cfg$model_type %in% c("sir_sus_frailty", "sir_trans_frailty")) {
    # Per-bin vac counts from the allocation_seed; matches build_frailty_mod
    # in estimate_from_data.R. Without this, run_stoch_frailty_cd would
    # default to round(f*n_total) per bin (no allocation variability).
    sd_pop  <- if (cfg$sd > 0) cfg$sd else cfg$sd_trans
    fr      <- get_frailty(sd = sd_pop, n = cfg$n_frailty)
    n_total <- round((cfg$N_cont + cfg$N_vac) * fr$p)
    bin     <- rep(seq_along(n_total), n_total)
    set.seed(cfg$allocation_seed)
    cfg$.vac_counts <- tabulate(bin[sample(length(bin), cfg$N_vac)],
                                nbins = length(n_total))
    set.seed(NULL)
  } else if (cfg$model_type == "sir_multisite") {
    # Per-site vaccinated counts from the allocation_seed (fixed for the fit;
    # the EATE draws fresh site allocations). site_icc sets the dispersion.
    N_site_vec <- multisite_site_sizes(cfg$N_cont + cfg$N_vac, cfg$n_sites)
    set.seed(cfg$allocation_seed)
    cfg$.vac_sites <- multisite_vac_counts(
      N_site_vec, cfg$N_vac,
      cfg$N_vac / (cfg$N_cont + cfg$N_vac), cfg$site_icc)
    set.seed(NULL)
  }
  cfg
}

# ---------------------------------------------------------------------------
# Fit (grid search + Nelder-Mead with restarts)
# ---------------------------------------------------------------------------

# Negative-log kernel (synthetic) likelihood at (beta, alpha), CRN-averaged
# over the fit seeds. L_hat(theta) = mean_i K_h(C_i - d) is the Monte-Carlo
# likelihood of the observed count pair; a Gaussian product kernel of
# bandwidth h (count units) smooths the exact "fraction equal to d" so the
# surface is differentiable. Matches the data to the MODE (unbiased under a
# bimodal / low-seed final-size distribution), unlike mean-matching.
.kernel_negloglik <- function(out, d1, d2, h) {
  if (is.null(out) || nrow(out) == 0) return(NA_real_)
  k <- exp(-0.5 * (((out$C1 - d1) / h)^2 + ((out$C2 - d2) / h)^2))
  -log(mean(k) / (2 * pi * h * h) + 1e-300)
}

make_loss <- function(simulator, cfg) {
  # Common Random Numbers: the SAME fixed seeds are reused at every probe,
  # so the loss is deterministic and smooth in (beta, alpha). Averaging
  # over the seeds reduces the seed-specific bias of the optimum.
  seeds  <- as.integer(cfg$opt_seed) + seq_len(max(1L, cfg$n_seed_opt)) - 1L
  method <- cfg$fit_method %||% "mean"
  h      <- cfg$kernel_h %||% 3
  function(log_par) {
    par <- exp(log_par); beta <- par[1]; alpha <- par[2]
    vals <- vapply(seeds, function(s) {
      out <- tryCatch(simulator(beta, alpha, cfg$n_sim_opt, seed = s),
                      error = function(e) NULL)
      if (method == "kernel") return(.kernel_negloglik(out, cfg$data_C1, cfg$data_C2, h))
      if (is.null(out) || nrow(out) == 0) return(NA_real_)
      (mean(out$C1) - cfg$data_C1)^2 + (mean(out$C2) - cfg$data_C2)^2
    }, numeric(1))
    if (all(is.na(vals))) return(1e9)
    mean(vals, na.rm = TRUE)
  }
}

grid_search_start <- function(loss_fn, cfg) {
  best <- list(loss = Inf, log_par = c(0, 0))
  for (lb in seq(log_beta_lo, log_beta_hi, length.out = cfg$grid_n)) {
    for (la in seq(log_alpha_lo, log_alpha_hi, length.out = cfg$grid_n)) {
      l <- loss_fn(c(lb, la))
      if (is.finite(l) && l < best$loss) best <- list(loss = l, log_par = c(lb, la))
    }
  }
  best
}

fit_one <- function(simulator, cfg) {
  loss  <- make_loss(simulator, cfg)
  start <- grid_search_start(loss, cfg)
  o <- optim(par = start$log_par, fn = loss, method = "Nelder-Mead",
             control = list(maxit = cfg$optim_maxit, reltol = 1e-4))
  best <- o
  for (rs in seq_len(cfg$n_restarts)) {
    if (o$convergence == 0) break
    if (best$value < cfg$restart_loss_threshold) break
    o <- optim(par = o$par, fn = loss, method = "Nelder-Mead",
               control = list(maxit = cfg$optim_maxit, reltol = 1e-4))
    if (o$value < best$value) best <- o
  }
  pars <- exp(best$par)
  list(beta = pars[1], alpha = pars[2],
       loss = best$value, convergence = best$convergence,
       grid_start_beta = exp(start$log_par[1]),
       grid_start_alpha = exp(start$log_par[2]),
       grid_start_loss = start$loss)
}

# Kernel-likelihood posterior covariance: the Laplace covariance is the inverse
# observed information = [Hessian of -log L_hat]^-1 at the MLE. The kernel
# broadens the data-space covariance by H = diag(h^2), inflating the posterior;
# rather than subtract it (unreliable under bimodality, since the mode-mass term
# dominates), extrapolate Var(h) -> 0 linearly in h^2 over `hs`, which keeps
# that information. Uses CRN (one fixed seed) and a larger n_sim, since the
# Hessian of a Monte-Carlo likelihood is noisy. Falls back to NA on failure.
kernel_posterior_cov <- function(simulator, cfg, beta, alpha) {
  d1 <- cfg$data_C1; d2 <- cfg$data_C2
  seed <- cfg$post_cov_seed
  nsim <- cfg$kernel_hess_nsim %||% 8000L
  hs   <- cfg$kernel_hess_hs   %||% c(2, 3, 4)
  nll <- function(b, a, h)
    .kernel_negloglik(simulator(b, a, nsim, seed = seed), d1, d2, h)
  h_rel <- cfg$kernel_hess_hrel %||% 0.03
  hess_cov_at <- function(h) {
    hb <- beta * h_rel; ha <- alpha * h_rel
    if (beta  - hb <= 0) hb <- beta  / 2
    if (alpha - ha <= 0) ha <- alpha / 2
    f0  <- nll(beta, alpha, h)
    fbb <- (nll(beta + hb, alpha, h) - 2 * f0 + nll(beta - hb, alpha, h)) / hb^2
    faa <- (nll(beta, alpha + ha, h) - 2 * f0 + nll(beta, alpha - ha, h)) / ha^2
    fba <- (nll(beta + hb, alpha + ha, h) - nll(beta + hb, alpha - ha, h) -
            nll(beta - hb, alpha + ha, h) + nll(beta - hb, alpha - ha, h)) / (4 * hb * ha)
    tryCatch(solve(matrix(c(fbb, fba, fba, faa), 2)),
             error = function(e) matrix(NA_real_, 2, 2))
  }
  covs <- lapply(hs, hess_cov_at)
  h2   <- hs^2
  extr <- function(idx) {
    y <- vapply(covs, function(C) C[idx], numeric(1))
    if (any(!is.finite(y))) return(NA_real_)
    unname(coef(stats::lm(y ~ h2))[1])                 # intercept = h -> 0
  }
  cov0 <- matrix(vapply(1:4, extr, numeric(1)), 2)
  cov0[1, 2] <- cov0[2, 1] <- 0.5 * (cov0[1, 2] + cov0[2, 1])
  dimnames(cov0) <- list(c("beta", "alpha"), c("beta", "alpha"))
  list(cov = cov0,
       sd  = sqrt(pmax(c(beta = cov0[1, 1], alpha = cov0[2, 2]), 0)))
}

# ---------------------------------------------------------------------------
# VE at the fitted (beta, alpha)
# ---------------------------------------------------------------------------

compute_ve <- function(cfg, beta, alpha) {
  N_total  <- cfg$N_cont + cfg$N_vac
  vac_frac <- cfg$N_vac / N_total
  sus      <- c(1, alpha)
  tp       <- seq(1, cfg$t_star, 1)

  ve <- switch(cfg$model_type,
    linear = get_stoch_eate_linear(
      beta = beta, susceptibility = sus, f = vac_frac, N = N_total,
      t = cfg$t_star, n_vac = cfg$ve_n_vac, n_rep = cfg$ve_n_rep,
      dt = cfg$dt, timepoints = tp, mc.cores = cfg$inner_cores),
    sir = get_stoch_eate_sir(
      beta = beta, susceptibility = sus, f = vac_frac, N = N_total,
      t = cfg$t_star, gamma = cfg$gamma, I_ini = cfg$I_ini_2g,
      n_vac = cfg$ve_n_vac, n_rep = cfg$ve_n_rep,
      dt = cfg$dt, timepoints = tp, mc.cores = cfg$inner_cores),
    sir_multisite = get_stoch_eate_sir_multisite(
      beta = beta, susceptibility = sus, f = vac_frac, N = N_total,
      t = cfg$t_star, gamma = cfg$gamma, I_ini = cfg$I_ini_2g,
      n_sites = cfg$n_sites, site_icc = cfg$site_icc,
      n_vac = cfg$ve_n_vac, n_rep = cfg$ve_n_rep,
      dt = cfg$dt, timepoints = tp, mc.cores = cfg$inner_cores),
    sir_sus_frailty = ,
    sir_trans_frailty = get_stoch_eate_frailty(
      alpha = alpha, sd = cfg$sd, sd_trans = cfg$sd_trans, beta = beta,
      f = vac_frac, N = N_total / 2, t = cfg$t_star,
      n_frailty = cfg$n_frailty, gamma = cfg$gamma,
      I_ini_total = sum(cfg$I_ini_2g),
      n_vac = cfg$ve_n_vac, n_rep = cfg$ve_n_rep,
      timepoints = tp, mc.cores = cfg$inner_cores,
      frailty_amp = cfg$frailty_amp %||% 2.5),
    network = get_stoch_eate_network(
      beta = beta, susceptibility = sus, f = vac_frac, N = N_total,
      t = cfg$t_star, c_ij = cfg$.c_ij, k_mean = cfg$mean_k,
      gamma = cfg$gamma, n_vac = cfg$ve_n_vac, n_rep = cfg$ve_n_rep,
      timepoints = tp, init_I = cfg$init_I_nw,
      mc.cores = cfg$inner_cores),
    stop("Unknown model_type: ", cfg$model_type))
  setDT(ve)
  ve[, model := cfg$model_type]
  ve[]
}

# ---------------------------------------------------------------------------
# Posterior sampling + VE with propagated parameter uncertainty
# ---------------------------------------------------------------------------

# Draw K samples (beta, alpha) from MVN(c(beta_hat, alpha_hat), cov), via
# Cholesky. Reject samples that fall into the negative orthant (rare when
# posterior is tight; if the rejection fraction is high it signals the
# Gaussian approximation is wrong for that fit). Returns a [<=K, 2] matrix.
sample_posterior <- function(beta_hat, alpha_hat, cov, K,
                             max_attempts_factor = 5) {
  mu <- c(beta_hat, alpha_hat)
  L  <- tryCatch(chol(cov),
                 error = function(e) chol(cov + diag(1e-10, 2)))
  out      <- matrix(NA_real_, nrow = K, ncol = 2,
                     dimnames = list(NULL, c("beta", "alpha")))
  n_filled <- 0L
  attempts <- 0L
  max_attempts <- as.integer(max_attempts_factor * K)
  while (n_filled < K && attempts < max_attempts) {
    chunk <- K - n_filled
    Z       <- matrix(rnorm(2L * chunk), nrow = chunk, ncol = 2)
    samples <- sweep(Z %*% L, 2, mu, "+")
    valid   <- samples[, 1] > 0 & samples[, 2] > 0
    n_valid <- sum(valid)
    if (n_valid > 0L) {
      take <- min(n_valid, K - n_filled)
      out[(n_filled + 1L):(n_filled + take), ] <-
        samples[valid, , drop = FALSE][seq_len(take), , drop = FALSE]
      n_filled <- n_filled + take
    }
    attempts <- attempts + chunk
  }
  if (n_filled < K) {
    warning(sprintf("Only %d / %d valid posterior samples after %d attempts (Gaussian posterior may be a poor fit here)",
                    n_filled, K, attempts))
  }
  out[seq_len(n_filled), , drop = FALSE]
}

# For each (beta_k, alpha_k) sample, call compute_ve which returns the
# EATE function's full per-allocation, per-time output (long format with
# `sim` = inner allocation/replicate index). Tag with param_sample and
# the corresponding (beta_k, alpha_k) so downstream can decompose
# variance into "parameter" vs "allocation" components.
compute_ve_with_uncertainty <- function(cfg, fit, posterior_cov,
                                        K, n_rep_override, K_cores = 1L) {
  samples <- sample_posterior(fit$beta, fit$alpha, posterior_cov, K)
  cfg_u              <- cfg
  cfg_u$ve_n_rep     <- n_rep_override
  # Inside the K loop each worker owns one core — one dust call at a
  # time, single-threaded. Prevents nested oversubscription with K_cores
  # already using the full CPU budget.
  cfg_u$inner_cores  <- 1L
  rbindlist(parallel::mclapply(seq_len(nrow(samples)), function(k) {
    ve_k <- compute_ve(cfg_u, samples[k, "beta"], samples[k, "alpha"])
    ve_k[, param_sample := k]
    ve_k[, beta_k       := samples[k, "beta"]]
    ve_k[, alpha_k      := samples[k, "alpha"]]
    ve_k
  }, mc.cores = K_cores, mc.preschedule = FALSE), fill = TRUE)
}

# ---------------------------------------------------------------------------
# Per-job runner
# ---------------------------------------------------------------------------

run_one_job <- function(cfg) {
  cfg       <- materialise_cfg(cfg)

  # --- Fitting phase: one dust call at a time. Give it all the cores.
  cfg_fit                <- cfg
  cfg_fit$inner_cores    <- cores_per_node
  simulator              <- build_simulator(cfg_fit)

  message(glue("[{cfg$name}] fitting..."))
  fit <- fit_one(simulator, cfg_fit)
  message(glue("[{cfg$name}] fit: beta = {round(fit$beta, 4)}  alpha = {round(fit$alpha, 4)}  ",
               "loss = {round(fit$loss, 3)}  conv = {fit$convergence}"))

  message(glue("[{cfg$name}] posterior cov..."))
  # Always compute the moment pieces (base, Sigma, J) for the misfit
  # diagnostics below. For the kernel fit the posterior cov/sd come from the
  # kernel-Hessian (h->0 extrapolated) instead of the moment sandwich, since
  # the fit matched the mode, not the mean.
  pcov <- estimate_posterior_cov(simulator, fit$beta, fit$alpha,
                                 n_sim = cfg$post_cov_n_sim,
                                 seed  = cfg$post_cov_seed,
                                 h_rel = cfg$post_cov_h)
  if ((cfg$fit_method %||% "mean") == "kernel") {
    kc <- kernel_posterior_cov(simulator, cfg_fit, fit$beta, fit$alpha)
    pcov$cov <- kc$cov
    pcov$sd  <- kc$sd
  }
  message(glue("[{cfg$name}] sd_beta = {signif(pcov$sd['beta'], 3)}  ",
               "sd_alpha = {signif(pcov$sd['alpha'], 3)}"))

  # Fit-quality diagnostics. Sigma is the per-realisation Cov(C); residuals
  # are the model-vs-data gap at the optimum (pcov$base, post_cov_n_sim reps).
  #   loss_floor = E[fit$loss] at a perfect fit = tr(Sigma) / n_sim_opt.
  #   loss_chisq = Mahalanobis distance of the data from the model mean using
  #                the PER-REALISATION Sigma: "is (data_C1,data_C2) a plausible
  #                draw from the model?" (~chi-square with 2 df; > ~6 = the
  #                model's reachable mean is implausibly far from the data =
  #                genuine misfit). NB: dividing by Sigma/n instead would test
  #                whether the mean matches the data to ~1/sqrt(n) of an SD,
  #                which over-flags even good fits and low-variance models.
  Sig        <- pcov$Sigma
  resid      <- c(mean(pcov$base$C1) - cfg$data_C1,
                  mean(pcov$base$C2) - cfg$data_C2)
  resid_C1   <- resid[1]
  resid_C2   <- resid[2]
  loss_floor <- (Sig[1, 1] + Sig[2, 2]) / cfg$n_sim_opt
  loss_chisq <- tryCatch(as.numeric(t(resid) %*% solve(Sig) %*% resid),
                         error = function(e) NA_real_)
  message(glue("[{cfg$name}] loss = {round(fit$loss, 3)}  ",
               "floor = {signif(loss_floor, 3)}  ",
               "chisq = {signif(loss_chisq, 3)}",
               "{if (is.finite(loss_chisq) && loss_chisq > 6) '  [MISFIT]' else ''}"))

  # Point-estimate VE (single call, one lot of n_vac x n_rep dust
  # replicates). Same phase strategy as fit.
  message(glue("[{cfg$name}] VE..."))
  ve <- compute_ve(cfg_fit, fit$beta, fit$alpha)

  # VE with uncertainty: K sequential compute_ve calls dominate the
  # total time. Run them in parallel across cores_per_node; each K
  # worker uses inner_cores = 1 (set inside the helper).
  message(glue("[{cfg$name}] VE with uncertainty (K = {cfg$K_post_samples}, ",
               "K_cores = {cores_per_node})..."))
  ve_unc <- compute_ve_with_uncertainty(cfg, fit, pcov$cov,
                                        K = cfg$K_post_samples,
                                        n_rep_override = cfg$ve_n_rep_uncert,
                                        K_cores = cores_per_node)

  list(
    name            = cfg$name,
    experiment_id   = cfg$experiment_id   %||% NA_character_,
    model_type      = cfg$model_type,
    fit_method      = cfg$fit_method      %||% "mean",
    pl_alpha        = cfg$pl_alpha        %||% NA_real_,
    network_seed    = cfg$network_seed    %||% NA_integer_,
    allocation_seed = cfg$allocation_seed %||% NA_integer_,
    fit             = fit,
    posterior_cov   = list(cov = pcov$cov, J = pcov$J,
                           Sigma = pcov$Sigma, sd = pcov$sd),
    loss_floor      = loss_floor,
    loss_chisq      = loss_chisq,
    resid_C1        = resid_C1,
    resid_C2        = resid_C2,
    ve              = ve,
    ve_uncertainty  = ve_unc
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Slurm array dispatch
# ---------------------------------------------------------------------------
# Shuffle configs once, split into `nodes` chunks, this slurm task
# processes its chunk. Configs inside a chunk run SEQUENTIALLY —
# parallelism is entirely inside each config (fit dust threading +
# K-loop parallelism), driven by `cores_per_node`. No outer mclapply
# to avoid nested-fork thrashing.

id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (id == "") id <- 1
id <- as.numeric(id)

chunk2  <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))
set.seed(1)
chunked <- chunk2(sample(seq_along(configs)), nodes)
set.seed(NULL)

ids <- as.integer(chunked[[id]])
message(glue("Task {id}/{nodes}: running {length(ids)} configs ",
             "({paste(ids, collapse=',')}), cores_per_node = {cores_per_node}"))

run_id <- function(i) {
  message(glue("Start {i}: {configs[[i]]$name}"))
  out <- tryCatch(run_one_job(configs[[i]]),
                  error = function(e) { message("Error: ", conditionMessage(e)); NULL })
  message(glue("End {i}"))
  out
}

results <- lapply(ids, run_id)

saveRDS(results, file.path(out_dir, glue("results_{id}.RDS")))
message(glue("Wrote {out_dir}/results_{id}.RDS"))
