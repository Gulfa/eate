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

# Sweep the homogeneous SIR model over multiple initial-infected settings
# WITHIN each experiment; each appears as its OWN model ("sir_i<total>") in the
# plots, so the I_ini effect is directly comparable. NULL = a single plain
# "sir" at the experiment's I_ini_2g. Each entry is a c(unvac, vac) vector.
sir_I_inis <- NULL          # e.g. list(c(1, 1), c(5, 5), c(10, 10))

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
kernel_h_frac    <- 0.1              # kernel bandwidth as a fraction of each data
                                     # count (scales with N; floored at 1 count)
kernel_hess_mult <- c(0.6, 0.8, 1.0, 1.2, 1.5) # bandwidth multipliers for ->0 extrapolation
kernel_hess_nsim <- 8000L            # sims for the (noisy) kernel Hessian
kernel_hess_hrel <- 0.03             # FD step (relative) for the 2nd-derivative Hessian

# Grid likelihood surface (the kernel posterior; replaces the FD Hessian).
# n x n nodes around the fit, box +/- span in log units (auto-expanded until
# the likelihood-ratio region is interior). Gives the MLE, profile LR
# intervals, moments, and the draws propagated into VE.
grid_post_n    <- 25L
grid_post_nsim <- 4000L
grid_post_span <- 0.7

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

# Two-block effect-modification model. split_frac = fraction of the population
# in compartment A; the vaccinated susceptibility is alpha in A and
# split_alpha_prod/alpha in B (so alpha_A*alpha_B = split_alpha_prod). One
# index case per realisation seeds a random compartment. split_init_I = index
# cases placed in the seeded compartment.
split_frac       <- 0.75
split_alpha_prod <- 0.5
split_init_I     <- 1

base_common <- list(
  gamma = gamma, dt = dt,
  n_sim_opt = n_sim_opt, optim_maxit = optim_maxit,
  n_restarts = n_restarts, restart_loss_threshold = restart_loss_threshold,
  grid_n = grid_n, opt_seed = opt_seed, n_seed_opt = n_seed_opt,
  # inner_cores is set by run_one_job per phase from cores_per_node.
  inner_cores = 1L,
  post_cov_n_sim = post_cov_n_sim, post_cov_seed = post_cov_seed,
  post_cov_h = post_cov_h,
  fit_method = fit_method, kernel_h_frac = kernel_h_frac,
  kernel_hess_mult = kernel_hess_mult, kernel_hess_nsim = kernel_hess_nsim,
  kernel_hess_hrel = kernel_hess_hrel,
  grid_post_n = grid_post_n, grid_post_nsim = grid_post_nsim,
  grid_post_span = grid_post_span,
  split_frac = split_frac, split_alpha_prod = split_alpha_prod,
  split_init_I = split_init_I,
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
  # Homogeneous SIR: either a single plain "sir", or one config per I_ini in
  # sir_I_inis, each shown as its own model "sir_i<total>" (sim_type stays
  # "sir" so the simulator/EATE switch is unchanged).
  sir_variants <- if (is.null(sir_I_inis)) list(NULL) else sir_I_inis
  for (iv in sir_variants) {
    if (is.null(iv)) {
      cs[[length(cs)+1]] <- modifyList(base, list(
        name = glue("{exp$id}__sir"),
        model_type = "sir", ve_n_vac = 1))
    } else {
      tag <- sprintf("i%02d", sum(iv))
      cs[[length(cs)+1]] <- modifyList(base, list(
        name = glue("{exp$id}__sir_{tag}"),
        model_type = glue("sir_{tag}"), sim_type = "sir",
        ve_n_vac = 1, I_ini_2g = iv))
    }
  }
  # Two-block effect-modification model: same vaccine, different effect in two
  # non-mixing compartments; one random index case decides which ignites, so
  # mode != mean in the vaccine ratio. Randomisation is simple within
  # compartments (exchangeable), so a single config + ve_n_vac = 1.
  # Fit with the KERNEL (mode) method: the moment-matching sandwich
  # J^-1 Sigma J^-T breaks for this model -- at large N the two compartments
  # ignite near-deterministically, so per-realisation (C1,C2) collapses to two
  # points (rank-1 Sigma) -> non-PD covariance / sd = 0; and alpha enters the
  # two compartments oppositely (alpha vs prod/alpha) so its column of J
  # nearly cancels -> wide/degenerate alpha. The kernel-Hessian posterior at
  # the ignited compartment's mode avoids both.
  cs[[length(cs)+1]] <- modifyList(base, list(
    name = glue("{exp$id}__sir_split_effect"),
    model_type = "sir_split_effect", ve_n_vac = 1,
    fit_method = "kernel"))
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
  # sim_type drives the simulator; model_type may be a display variant
  # (e.g. "sir_i02") that shares a simulator.
  switch(cfg$sim_type %||% cfg$model_type,
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
    sir_split_effect = function(beta, alpha, n_sim, seed = NULL) {
      # Two non-mixing compartments, vaccinated susceptibility alpha in A and
      # split_alpha_prod/alpha in B. One index case per realisation lands in a
      # random compartment (stratified: split_frac of particles seeded in A,
      # rest in B). Pool arms across compartments: C1 = unvac cases, C2 = vac.
      f  <- cfg$N_vac / (cfg$N_cont + cfg$N_vac)
      sp <- .split_effect_setup(alpha, f, cfg$N_cont + cfg$N_vac,
                                cfg$split_frac, cfg$split_alpha_prod)
      I0 <- as.integer(max(1, round(cfg$split_init_I %||% 1)))
      nA <- round(cfg$split_frac * n_sim)
      run_b <- function(gseed, nsim, sd) {
        if (nsim <= 0) return(NULL)
        Ii <- integer(4); Ii[gseed] <- I0
        o <- run_stoch_cd_dust(sp$mm, beta = beta, N = sp$Ngrp, t = cfg$t_star,
                               I_ini = Ii, susceptibility = sp$sus,
                               gamma = cfg$gamma, dt = cfg$dt,
                               timepoints = seq(1, cfg$t_star, 1),
                               n_sim = nsim, cores = cfg$inner_cores, seed = sd)
        setDT(o); fin <- o[time == cfg$t_star]
        data.table(C1 = fin$C1 + fin$C3, C2 = fin$C2 + fin$C4)
      }
      sdB <- if (is.null(seed)) NULL else as.integer(seed) + 1L
      rbindlist(list(run_b(1L, nA, seed), run_b(3L, n_sim - nA, sdB)))
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
# likelihood of the observed count pair; a Gaussian product kernel of per-arm
# bandwidth (h1, h2) smooths the exact "fraction equal to d" so the surface is
# differentiable. Matches the data to the MODE (unbiased under a bimodal /
# low-seed final-size distribution), unlike mean-matching.
.kernel_negloglik <- function(out, d1, d2, h1, h2) {
  if (is.null(out) || nrow(out) == 0) return(NA_real_)
  k <- exp(-0.5 * (((out$C1 - d1) / h1)^2 + ((out$C2 - d2) / h2)^2))
  -log(mean(k) / (2 * pi * h1 * h2) + 1e-300)
}

# Per-arm kernel bandwidths, scaled to the data magnitude so the kernel works
# at any N: h_arm = kernel_h_frac * data_arm (floored at 1 count). A fixed
# absolute bandwidth would be far too small once counts reach the thousands,
# leaving no simulated point within h of the data -> flat (L_hat=0) plateau.
.kernel_bw <- function(cfg) {
  frac <- cfg$kernel_h_frac %||% 0.1
  c(max(frac * cfg$data_C1, 1), max(frac * cfg$data_C2, 1))
}

make_loss <- function(simulator, cfg) {
  # Common Random Numbers: the SAME fixed seeds are reused at every probe,
  # so the loss is deterministic and smooth in (beta, alpha). Averaging
  # over the seeds reduces the seed-specific bias of the optimum.
  seeds  <- as.integer(cfg$opt_seed) + seq_len(max(1L, cfg$n_seed_opt)) - 1L
  method <- cfg$fit_method %||% "mean"
  bw     <- .kernel_bw(cfg)
  function(log_par) {
    par <- exp(log_par); beta <- par[1]; alpha <- par[2]
    vals <- vapply(seeds, function(s) {
      out <- tryCatch(simulator(beta, alpha, cfg$n_sim_opt, seed = s),
                      error = function(e) NULL)
      if (method == "kernel") return(.kernel_negloglik(out, cfg$data_C1, cfg$data_C2, bw[1], bw[2]))
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

# ---------------------------------------------------------------------------
# Grid likelihood surface (replaces the finite-difference Hessian)
# ---------------------------------------------------------------------------
# Evaluate -log L_hat on a grid around the fit under CRN. With only 2
# parameters this is cheap and fully parallel, and -- unlike a numerical
# Hessian -- it does NOT differentiate a Monte-Carlo function: FD amplifies MC
# noise by ~1/h^2 (which produced sd=0 / non-PD / unstable rho, with the
# across-seed CV of sd_beta measured at 0.4-1.2), whereas summing over a grid
# averages it. From one surface we get:
#   - the MLE (argmax) as a check on the optimiser,
#   - frequentist likelihood-ratio regions: {2*(logLmax - logL) <= chisq},
#     profiled per parameter (chisq_1 = 3.84 at 95%),
#   - moments and DRAWS by weighting nodes by L (flat prior over the box),
#     used to propagate parameter uncertainty into VE.
# Box is set in log space around the fit and expanded until the region of
# interest is interior; nodes are spaced linearly so a flat prior in (beta,
# alpha) makes the weights proportional to L.
#
# NB the chi-square calibration is exact for a true likelihood; L_hat is
# kernel-smoothed (bandwidth h) and MC-estimated, so intervals are h-inflated
# and only approximately calibrated.
grid_posterior <- function(simulator, cfg, beta_hat, alpha_hat, cores = 1L) {
  n_grid <- cfg$grid_post_n    %||% 25L
  nsim   <- cfg$grid_post_nsim %||% 4000L
  span   <- cfg$grid_post_span %||% 0.7      # +/- in log units, then expanded
  seed   <- cfg$post_cov_seed
  bw     <- .kernel_bw(cfg)
  d1     <- cfg$data_C1; d2 <- cfg$data_C2
  cut_in <- 12                                # region that must be interior

  eval_grid <- function(bs, as_) {
    g <- as.data.table(expand.grid(beta = bs, alpha = as_))
    v <- unlist(parallel::mclapply(seq_len(nrow(g)), function(i)
      .kernel_negloglik(simulator(g$beta[i], g$alpha[i], nsim, seed = seed),
                        d1, d2, bw[1], bw[2]),
      mc.cores = cores, mc.preschedule = TRUE))
    g[, nll := as.numeric(v)][]
  }

  lb <- log(beta_hat); la <- log(alpha_hat)
  g <- NULL
  for (it in seq_len(3L)) {
    bs  <- seq(exp(lb - span), exp(lb + span), length.out = n_grid)
    as_ <- seq(exp(la - span), exp(la + span), length.out = n_grid)
    g   <- eval_grid(bs, as_)
    if (!any(is.finite(g$nll))) return(NULL)
    g[, dev := 2 * (nll - min(nll, na.rm = TRUE))]
    edge <- g[beta %in% range(bs) | alpha %in% range(as_)]
    if (!any(is.finite(edge$dev)) || min(edge$dev, na.rm = TRUE) > cut_in) break
    span <- span * 1.6                        # region hits the edge -> widen
  }
  g <- g[is.finite(nll)]
  if (!nrow(g)) return(NULL)
  g[, dev := 2 * (nll - min(nll))]
  g[, w := exp(-(nll - min(nll)))]
  g[, w := w / sum(w)]

  mu_b <- sum(g$w * g$beta); mu_a <- sum(g$w * g$alpha)
  vb   <- sum(g$w * (g$beta  - mu_b)^2)
  va   <- sum(g$w * (g$alpha - mu_a)^2)
  cba  <- sum(g$w * (g$beta - mu_b) * (g$alpha - mu_a))
  cov0 <- matrix(c(vb, cba, cba, va), 2,
                 dimnames = list(c("beta", "alpha"), c("beta", "alpha")))

  # Profile likelihood-ratio intervals (chisq_1 = 3.84 at 95%).
  prof <- function(par) {
    p <- g[, .(dev = min(dev)), by = par][dev <= 3.84]
    if (!nrow(p)) return(c(NA_real_, NA_real_))
    range(p[[par]])
  }
  mle <- g[which.min(nll)]
  list(cov = cov0,
       sd  = sqrt(c(beta = vb, alpha = va)),
       mean = c(beta = mu_b, alpha = mu_a),
       mle  = c(beta = mle$beta, alpha = mle$alpha),
       ci_beta = prof("beta"), ci_alpha = prof("alpha"),
       span = span, n_grid = n_grid,
       d_beta = diff(sort(unique(g$beta))[1:2]),
       d_alpha = diff(sort(unique(g$alpha))[1:2]),
       grid = g[, .(beta, alpha, nll, dev, w)])
}

# Draw K parameter pairs by weighting grid nodes by L (flat prior over the
# box), jittered within a cell so the draws are not stuck on grid lines.
sample_grid_posterior <- function(gp, K) {
  if (is.null(gp) || !nrow(gp$grid)) return(NULL)
  g   <- gp$grid
  idx <- sample.int(nrow(g), K, replace = TRUE, prob = g$w)
  out <- cbind(beta  = g$beta[idx]  + runif(K, -gp$d_beta  / 2, gp$d_beta  / 2),
               alpha = g$alpha[idx] + runif(K, -gp$d_alpha / 2, gp$d_alpha / 2))
  out[, 1] <- pmax(out[, 1], 1e-8)
  out[, 2] <- pmax(out[, 2], 1e-8)
  colnames(out) <- c("beta", "alpha")
  out
}

# Kernel-likelihood posterior covariance: the Laplace covariance is the inverse
# observed information = [Hessian of -log L_hat]^-1 at the MLE. The kernel
# broadens the data-space covariance by H = diag(h^2), inflating the posterior;
# rather than subtract it (unreliable under bimodality, since the mode-mass term
# dominates), extrapolate the bandwidth -> 0 (variance linear in the squared
# scale multiplier), which keeps that information. The per-arm base bandwidth
# is scaled by the multipliers in kernel_hess_mult. Uses CRN (one fixed seed)
# and a larger n_sim, since the Hessian of a MC likelihood is noisy.
kernel_posterior_cov <- function(simulator, cfg, beta, alpha) {
  d1 <- cfg$data_C1; d2 <- cfg$data_C2
  seed <- cfg$post_cov_seed
  nsim <- cfg$kernel_hess_nsim %||% 8000L
  bw   <- .kernel_bw(cfg)                               # per-arm base bandwidth
  mult <- cfg$kernel_hess_mult %||% c(0.7, 1.0, 1.4)    # scale factors
  nll <- function(b, a, m)
    .kernel_negloglik(simulator(b, a, nsim, seed = seed), d1, d2, m * bw[1], m * bw[2])
  h_rel     <- cfg$kernel_hess_hrel %||% 0.03
  floor_nll <- -log(1e-300)                             # value when L_hat = 0
  # Adaptive FD step: shrink while a +/- perturbation leaves the kernel (nll
  # near the floor, or a huge jump) -> keeps the 2nd difference in the
  # informative region; grow while the perturbation is too small to register.
  # A fixed relative step overshoots the SHARP likelihood at large N (a small
  # param change moves the near-deterministic outbreak past the bandwidth),
  # hitting the floor -> spuriously huge Hessian -> sd = 0.
  find_step <- function(f_pert, base_step, f0) {
    h <- base_step
    for (it in 1:16) {
      fp <- f_pert(h); fm <- f_pert(-h)
      if (!is.finite(fp) || !is.finite(fm)) { h <- h / 2; next }
      dmax <- max(fp, fm) - f0
      if (max(fp, fm) > floor_nll - 5 || dmax > 15) { h <- h / 2; next }  # overshoot
      if (dmax < 0.05) { h <- h * 2; next }                              # understep
      break
    }
    h
  }
  # Hessian of -log L_hat at (beta, alpha) for bandwidth multiplier m, by
  # fitting a local QUADRATIC surface via least squares over a small design
  # that is kept INSIDE the kernel. This estimates the cross term (and hence
  # the beta-alpha correlation) stably -- the old corner second-difference
  # moved both params at once, so its corners left the kernel and the cross
  # term was garbage (|rho| >= 1 -> non-PD). Design radii come from the
  # per-axis adaptive steps, shrunk together until the corners are in-kernel.
  hessian_at <- function(m) {
    f0 <- nll(beta, alpha, m)
    if (!is.finite(f0) || f0 > floor_nll - 5) return(NULL)
    hb <- find_step(function(dh) nll(beta + dh, alpha, m), beta  * h_rel, f0)
    ha <- find_step(function(dh) nll(beta, alpha + dh, m), alpha * h_rel, f0)
    hb <- min(hb, beta / 2); ha <- min(ha, alpha / 2)
    dirs <- rbind(c(1,0), c(-1,0), c(0,1), c(0,-1),
                  c(1,1), c(-1,-1), c(1,-1), c(-1,1))
    for (shrink in 1:6) {                              # keep corners in-kernel
      pts <- cbind(beta + dirs[,1]*hb, alpha + dirs[,2]*ha)
      y   <- vapply(seq_len(nrow(pts)),
                    function(i) nll(pts[i,1], pts[i,2], m), numeric(1))
      if (all(is.finite(y)) && max(y) < floor_nll - 5) break
      hb <- hb * 0.7; ha <- ha * 0.7
    }
    if (!all(is.finite(y))) return(NULL)
    db <- pts[,1] - beta; da <- pts[,2] - alpha
    X  <- cbind(db, da, db^2, da^2, db*da)             # centre absorbs f0
    cf <- tryCatch(unname(coef(lm(I(y - f0) ~ X - 1))), error = function(e) NULL)
    if (is.null(cf) || any(!is.finite(cf))) return(NULL)
    matrix(c(2*cf[3], cf[5], cf[5], 2*cf[4]), 2)       # H = curvature
  }
  hess_cov_at <- function(m) {
    H <- hessian_at(m)
    if (is.null(H)) return(matrix(NA_real_, 2, 2))
    tryCatch(solve(H), error = function(e) matrix(NA_real_, 2, 2))
  }
  covs <- lapply(mult, hess_cov_at)
  m2   <- mult^2
  imin <- which.min(mult)                              # least bandwidth inflation
  extr <- function(idx) {
    y <- vapply(covs, function(C) C[idx], numeric(1))
    if (any(!is.finite(y))) return(covs[[imin]][idx])  # fall back to smallest bw
    unname(coef(stats::lm(y ~ m2))[1])                 # intercept = bandwidth -> 0
  }
  cov0 <- matrix(vapply(1:4, extr, numeric(1)), 2)
  cov0[1, 2] <- cov0[2, 1] <- 0.5 * (cov0[1, 2] + cov0[2, 1])
  # Project to the nearest positive-definite matrix (eigenvalue floor). Handles
  # a residual noisy/ridge cross term or a variance that extrapolated <= 0,
  # without pretending a specific correlation.
  if (any(!is.finite(cov0))) cov0 <- covs[[imin]]
  if (any(!is.finite(cov0))) cov0 <- diag(c(1e-4, 1e-4))
  ev  <- eigen(cov0, symmetric = TRUE)
  lam <- pmax(ev$values, max(ev$values, 0) * 1e-3 + 1e-12)
  cov0 <- ev$vectors %*% diag(lam) %*% t(ev$vectors)
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

  ve <- switch(cfg$sim_type %||% cfg$model_type,
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
    sir_split_effect = get_stoch_eate_sir_split_effect(
      beta = beta, susceptibility = sus, f = vac_frac, N = N_total,
      t = cfg$t_star, gamma = cfg$gamma, I_ini_total = cfg$split_init_I %||% 1,
      split_frac = cfg$split_frac, split_alpha_prod = cfg$split_alpha_prod,
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
  mu    <- c(beta_hat, alpha_hat)
  empty <- matrix(numeric(0), nrow = 0, ncol = 2,
                  dimnames = list(NULL, c("beta", "alpha")))
  # Guard: a singular Jacobian (e.g. a fit pinned at a parameter bound)
  # yields an NA / non-positive-definite covariance. Return no samples so the
  # job degrades gracefully (empty VE-with-uncertainty) instead of crashing.
  if (any(!is.finite(cov)) || any(!is.finite(mu))) {
    warning("Non-finite posterior covariance/mean; skipping VE-with-uncertainty.")
    return(empty)
  }
  L <- tryCatch(chol(cov),
                error = function(e) tryCatch(chol(cov + diag(1e-10, 2)),
                                             error = function(e2) NULL))
  if (is.null(L) || any(!is.finite(L))) {
    warning("Posterior covariance not positive-definite; skipping VE-with-uncertainty.")
    return(empty)
  }
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
                                        K, n_rep_override, K_cores = 1L,
                                        samples = NULL) {
  # `samples` (e.g. grid-posterior draws) takes precedence; otherwise fall
  # back to the Gaussian MVN(fit, posterior_cov) approximation.
  if (is.null(samples))
    samples <- sample_posterior(fit$beta, fit$alpha, posterior_cov, K)
  if (is.null(samples) || !nrow(samples)) return(data.table())
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
  # A fit pinned at a parameter bound usually means the data target is
  # unreachable for this setup (e.g. it sits in the bimodal fizzle/takeoff gap).
  .at_bnd <- function(x, lo, hi) {
    span <- abs(hi - lo)
    abs(log(x) - lo) < 0.02 * span || abs(log(x) - hi) < 0.02 * span
  }
  if (.at_bnd(fit$beta, log_beta_lo, log_beta_hi) ||
      .at_bnd(fit$alpha, log_alpha_lo, log_alpha_hi))
    message(glue("[{cfg$name}] WARNING: fit pinned at a parameter bound -- ",
                 "target likely infeasible for this N / seed setup."))

  message(glue("[{cfg$name}] posterior cov..."))
  # Always compute the moment pieces (base, Sigma, J) for the misfit
  # diagnostics below. For the kernel fit the posterior comes from the GRID
  # likelihood surface (no differentiation of a Monte-Carlo function), which
  # also supplies the draws used to propagate uncertainty into VE.
  pcov <- estimate_posterior_cov(simulator, fit$beta, fit$alpha,
                                 n_sim = cfg$post_cov_n_sim,
                                 seed  = cfg$post_cov_seed,
                                 h_rel = cfg$post_cov_h)
  gpost <- NULL; post_draws <- NULL
  if ((cfg$fit_method %||% "mean") == "kernel") {
    cfg_g <- cfg_fit; cfg_g$inner_cores <- 1L      # parallelise over grid nodes
    gpost <- grid_posterior(build_simulator(cfg_g), cfg_g,
                            fit$beta, fit$alpha, cores = cores_per_node)
    if (!is.null(gpost)) {
      pcov$cov   <- gpost$cov
      pcov$sd    <- gpost$sd
      post_draws <- sample_grid_posterior(gpost, cfg$K_post_samples)
      message(glue("[{cfg$name}] grid MLE: beta = {signif(gpost$mle['beta'], 4)} ",
                   "alpha = {signif(gpost$mle['alpha'], 4)}  |  ",
                   "95% LR CI beta [{signif(gpost$ci_beta[1],3)}, {signif(gpost$ci_beta[2],3)}] ",
                   "alpha [{signif(gpost$ci_alpha[1],3)}, {signif(gpost$ci_alpha[2],3)}]"))
    } else {
      message(glue("[{cfg$name}] grid posterior failed; keeping moment sandwich."))
    }
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
                                        K_cores = cores_per_node,
                                        samples = post_draws)

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
    grid_post       = if (is.null(gpost)) NULL else
                      list(mle = gpost$mle, mean = gpost$mean,
                           ci_beta = gpost$ci_beta, ci_alpha = gpost$ci_alpha,
                           span = gpost$span, n_grid = gpost$n_grid),
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
#
# Set EATE_SOURCE_ONLY=1 to load the definitions (build_simulator,
# materialise_cfg, configs, ...) without running the array -- used by
# helper scripts such as compare_mode_mean.R.
if (Sys.getenv("EATE_SOURCE_ONLY") != "1") {

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

}   # end if (EATE_SOURCE_ONLY != "1")
