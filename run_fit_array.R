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
# Shared targets / fit / VE knobs
# ---------------------------------------------------------------------------

N_cont    <- 100
N_vac     <- 100
data_C1   <- 50
data_C2   <- 25
t_star    <- 8

gamma     <- 1
I_ini_2g  <- c(10, 10)
init_I_nw <- 2
dt        <- 0.01

# Optim
n_sim_opt   <- 1000
optim_maxit <- 250
n_restarts  <- 3
restart_loss_threshold <- 5
grid_n      <- 6
log_beta_lo  <- log(0.01); log_beta_hi  <- log(5)
log_alpha_lo <- log(0.01); log_alpha_hi <- log(2)

# Posterior cov
post_cov_n_sim <- 1000
post_cov_seed  <- 1234L
post_cov_h     <- 0.01

# VE (stochastic EATE)
ve_n_vac <- 5
ve_n_rep <- 200

# Internal dust threading (per simulator call). Outer parallelism is via
# slurm array tasks, so keep this modest.
inner_cores <- 4

out_dir <- "output/fit_array_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Build the config grid:
#   - one config each for linear, sir, sir_sus_frailty, sir_trans_frailty
#   - n_networks x n_allocations configs for network
# ---------------------------------------------------------------------------

n_networks    <- 10
n_allocations <- 10
pl_alpha      <- 3
mean_k        <- 6

base <- list(
  N_cont = N_cont, N_vac = N_vac,
  data_C1 = data_C1, data_C2 = data_C2,
  t_star = t_star, gamma = gamma,
  I_ini_2g = I_ini_2g, init_I_nw = init_I_nw, dt = dt,
  n_sim_opt = n_sim_opt, optim_maxit = optim_maxit,
  n_restarts = n_restarts, restart_loss_threshold = restart_loss_threshold,
  grid_n = grid_n, inner_cores = inner_cores,
  post_cov_n_sim = post_cov_n_sim, post_cov_seed = post_cov_seed,
  post_cov_h = post_cov_h,
  ve_n_vac = ve_n_vac, ve_n_rep = ve_n_rep
)

configs <- list()

configs[[length(configs)+1]] <- modifyList(base, list(
  name = "linear", model_type = "linear"))

configs[[length(configs)+1]] <- modifyList(base, list(
  name = "sir", model_type = "sir"))

configs[[length(configs)+1]] <- modifyList(base, list(
  name = "sir_sus_frailty", model_type = "sir_sus_frailty",
  sd = 0.3, sd_trans = 0, n_frailty = 10))

configs[[length(configs)+1]] <- modifyList(base, list(
  name = "sir_trans_frailty", model_type = "sir_trans_frailty",
  sd = 0, sd_trans = 0.3, n_frailty = 10))

for (network_seed in seq_len(n_networks)) {
  for (alloc_seed in seq_len(n_allocations)) {
    configs[[length(configs)+1]] <- modifyList(base, list(
      name            = glue("network_n{network_seed}_a{alloc_seed}"),
      model_type      = "network",
      pl_alpha        = pl_alpha, mean_k = mean_k,
      network_seed    = network_seed,
      allocation_seed = alloc_seed))
  }
}

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
    sir_sus_frailty = function(beta, alpha, n_sim, seed = NULL) {
      at_tstar(run_stoch_frailty_cd(
        sd = cfg$sd, sd_trans = cfg$sd_trans,
        susceptibility = c(1, alpha), beta = beta,
        N = N_total / 2, t = cfg$t_star,
        n_frailty = cfg$n_frailty, gamma = cfg$gamma,
        I_ini_total = sum(cfg$I_ini_2g),
        timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores,
        method = "dust", dt = cfg$dt, f = 0.5, seed = seed), cfg$t_star)
    },
    sir_trans_frailty = function(beta, alpha, n_sim, seed = NULL) {
      at_tstar(run_stoch_frailty_cd(
        sd = cfg$sd, sd_trans = cfg$sd_trans,
        susceptibility = c(1, alpha), beta = beta,
        N = N_total / 2, t = cfg$t_star,
        n_frailty = cfg$n_frailty, gamma = cfg$gamma,
        I_ini_total = sum(cfg$I_ini_2g),
        timepoints = seq(1, cfg$t_star, 1),
        n_sim = n_sim, cores = cfg$inner_cores,
        method = "dust", dt = cfg$dt, f = 0.5, seed = seed), cfg$t_star)
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
  }
  cfg
}

# ---------------------------------------------------------------------------
# Fit (grid search + Nelder-Mead with restarts)
# ---------------------------------------------------------------------------

make_loss <- function(simulator, cfg) {
  function(log_par) {
    par <- exp(log_par); beta <- par[1]; alpha <- par[2]
    out <- tryCatch(simulator(beta, alpha, cfg$n_sim_opt),
                    error = function(e) NULL)
    if (is.null(out) || nrow(out) == 0) return(1e9)
    (mean(out$C1) - cfg$data_C1)^2 + (mean(out$C2) - cfg$data_C2)^2
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
    sir_sus_frailty = ,
    sir_trans_frailty = get_stoch_eate_frailty(
      alpha = alpha, sd = cfg$sd, sd_trans = cfg$sd_trans, beta = beta,
      f = vac_frac, N = N_total / 2, t = cfg$t_star,
      n_frailty = cfg$n_frailty, gamma = cfg$gamma,
      I_ini_total = sum(cfg$I_ini_2g),
      n_vac = cfg$ve_n_vac, n_rep = cfg$ve_n_rep,
      timepoints = tp, mc.cores = cfg$inner_cores),
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
# Per-job runner
# ---------------------------------------------------------------------------

run_one_job <- function(cfg) {
  cfg       <- materialise_cfg(cfg)
  simulator <- build_simulator(cfg)

  message(glue("[{cfg$name}] fitting..."))
  fit <- fit_one(simulator, cfg)
  message(glue("[{cfg$name}] fit: beta = {round(fit$beta, 4)}  alpha = {round(fit$alpha, 4)}  ",
               "loss = {round(fit$loss, 3)}  conv = {fit$convergence}"))

  message(glue("[{cfg$name}] posterior cov..."))
  pcov <- estimate_posterior_cov(simulator, fit$beta, fit$alpha,
                                 n_sim = cfg$post_cov_n_sim,
                                 seed  = cfg$post_cov_seed,
                                 h_rel = cfg$post_cov_h)
  message(glue("[{cfg$name}] sd_beta = {signif(pcov$sd['beta'], 3)}  ",
               "sd_alpha = {signif(pcov$sd['alpha'], 3)}"))

  message(glue("[{cfg$name}] VE..."))
  ve <- compute_ve(cfg, fit$beta, fit$alpha)

  list(
    name            = cfg$name,
    model_type      = cfg$model_type,
    network_seed    = cfg$network_seed    %||% NA_integer_,
    allocation_seed = cfg$allocation_seed %||% NA_integer_,
    fit             = fit,
    posterior_cov   = list(cov = pcov$cov, J = pcov$J,
                           Sigma = pcov$Sigma, sd = pcov$sd),
    ve              = ve
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Slurm array dispatch — same chunking pattern as run_array.R
# ---------------------------------------------------------------------------

N_nodes        <- 20
cores_per_node <- 4   # outer mclapply across jobs in this task's chunk

id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (id == "") id <- 1
id <- as.numeric(id)

chunk2  <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))
set.seed(1)
chunked <- chunk2(sample(seq_along(configs)), N_nodes)
set.seed(NULL)

ids <- as.integer(chunked[[id]])
message(glue("Task {id}: running {length(ids)} configs ({paste(ids, collapse=',')})"))

run_id <- function(i) {
  message(glue("Start {i}: {configs[[i]]$name}"))
  out <- tryCatch(run_one_job(configs[[i]]),
                  error = function(e) { message("Error: ", conditionMessage(e)); NULL })
  message(glue("End {i}"))
  out
}

results <- parallel::mclapply(ids, run_id,
                              mc.cores       = cores_per_node,
                              mc.preschedule = FALSE)

saveRDS(results, file.path(out_dir, glue("results_{id}.RDS")))
message(glue("Wrote {out_dir}/results_{id}.RDS"))
