# Stochastic VE sweep across models and vaccine-susceptibility values.
# Runs each model type at a fixed transmission rate (beta) for every alpha
# in `alphas`, keeping model-specific structural parameters as separate
# axes:
#   - linear / sir: no extra axis
#   - sir_sus_frailty, sir_trans_frailty: sweep over `frailty_sds`
#   - network: sweep over `pl_alphas` (Pareto exponent). n_network_seeds
#     independent contact graphs are averaged over per pl_alpha.
#
# For frailty and network the EATE function internally averages over
# `n_vac_allocs` vac allocations (per replicate) so no extra allocation
# loop is needed. n_rep dust replicates per allocation control the
# stochastic noise floor.
#
# Output: per-model VE(t) trajectories with 95% quantile bands, faceted
# by alpha in one plot and by model in another.
#
# Usage:
#   Rscript network_alpha_sweep.R

library(dplyr)
library(data.table)
library(ggplot2)
library(glue)

source("det_model.R")
source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Parameters (edit here)
# ---------------------------------------------------------------------------

# Vaccine-susceptibility values to sweep (each becomes one facet in the
# per-alpha plot).
alphas       <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# Network Pareto exponents (kept separate; each is its own curve in the
# network facet).
pl_alphas    <- c(2, 3, 5)

# Frailty SD values applied to both sus- and trans-frailty models. Each
# SD becomes its own curve.
frailty_sds  <- c(0.1, 0.3, 0.5)

# Common epidemic parameters
beta         <- 1.5   # transmission rate; R0 = beta/gamma in homog limit
gamma        <- 1
N            <- 200
t_horizon    <- 8
mean_k       <- 6
n_frailty    <- 10
init_I_2g    <- c(10, 10)   # initial infected per group for SIR/frailty
init_I_nw    <- 2

# EATE knobs
n_network_seeds <- 5    # independent Pareto graphs per pl_alpha
n_vac_allocs    <- 10   # inner allocation iterations (EATE n_vac)
n_rep           <- 100  # dust replicates per allocation

# Parallelism: outer mclapply across (model x alpha x variant) jobs;
# inner dust threading per call. Default sends all available cores to
# the outer loop (~115 jobs fit in 1-2 mclapply rounds on a 100-core
# machine, no slurm needed). Respects SLURM_CPUS_PER_TASK if set so a
# slurm sbatch with --cpus-per-task=N does the right thing.
outer_cores <- {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (nzchar(slurm_cpus)) as.integer(slurm_cpus)
  else max(1L, parallel::detectCores() - 1L)
}
inner_cores <- 1

out_dir <- "output/network_alpha_sweep"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

timepoints <- seq(1, t_horizon, 1)

# ---------------------------------------------------------------------------
# Materialise networks once per (pl_alpha, seed)
# ---------------------------------------------------------------------------

networks <- list()
for (pa in pl_alphas) {
  for (ns in seq_len(n_network_seeds)) {
    set.seed(1000L * round(pa * 10) + ns)
    key <- sprintf("pa%s_n%02d", format(pa, trim = TRUE), ns)
    networks[[key]] <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
  }
}
set.seed(NULL)
message(glue("Built {length(networks)} contact matrices ({length(pl_alphas)} pl_alphas x {n_network_seeds} seeds)"))

# ---------------------------------------------------------------------------
# Build the job list
# ---------------------------------------------------------------------------

jobs <- list()

for (a in alphas) {
  jobs[[length(jobs) + 1L]] <- list(model = "linear",  variant = "",
                                    alpha = a, sd_val = NA_real_)
  jobs[[length(jobs) + 1L]] <- list(model = "sir",     variant = "",
                                    alpha = a, sd_val = NA_real_)
  for (sd in frailty_sds) {
    jobs[[length(jobs) + 1L]] <- list(
      model = "sir_sus_frailty",
      variant = sprintf("sd%.2f", sd),
      alpha = a, sd_val = sd)
    jobs[[length(jobs) + 1L]] <- list(
      model = "sir_trans_frailty",
      variant = sprintf("sd%.2f", sd),
      alpha = a, sd_val = sd)
  }
  for (pa in pl_alphas) {
    for (ns in seq_len(n_network_seeds)) {
      jobs[[length(jobs) + 1L]] <- list(
        model = "network",
        variant = sprintf("pa%s", format(pa, trim = TRUE)),
        alpha = a, sd_val = NA_real_,
        pl_alpha = pa, network_seed = ns,
        key = sprintf("pa%s_n%02d", format(pa, trim = TRUE), ns))
    }
  }
}

message(glue("Built {length(jobs)} jobs total."))

# ---------------------------------------------------------------------------
# Per-job runner
# ---------------------------------------------------------------------------

run_job <- function(job) {
  susc <- c(1, job$alpha)
  ve <- switch(job$model,
    linear = get_stoch_eate_linear(
      beta = beta, susceptibility = susc, f = 0.5, N = N,
      t = t_horizon, n_vac = 1, n_rep = n_rep,
      timepoints = timepoints, mc.cores = inner_cores),
    sir = get_stoch_eate_sir(
      beta = beta, susceptibility = susc, f = 0.5, N = N,
      t = t_horizon, gamma = gamma, I_ini = init_I_2g,
      n_vac = 1, n_rep = n_rep,
      timepoints = timepoints, mc.cores = inner_cores),
    sir_sus_frailty = get_stoch_eate_frailty(
      alpha = job$alpha, sd = job$sd_val, sd_trans = 0, beta = beta,
      f = 0.5, N = N / 2, t = t_horizon,
      n_frailty = n_frailty, gamma = gamma,
      I_ini_total = sum(init_I_2g),
      n_vac = n_vac_allocs, n_rep = n_rep,
      timepoints = timepoints, mc.cores = inner_cores),
    sir_trans_frailty = get_stoch_eate_frailty(
      alpha = job$alpha, sd = 0, sd_trans = job$sd_val, beta = beta,
      f = 0.5, N = N / 2, t = t_horizon,
      n_frailty = n_frailty, gamma = gamma,
      I_ini_total = sum(init_I_2g),
      n_vac = n_vac_allocs, n_rep = n_rep,
      timepoints = timepoints, mc.cores = inner_cores),
    network = get_stoch_eate_network(
      beta = beta, susceptibility = susc, f = 0.5, N = N,
      t = t_horizon, c_ij = networks[[job$key]], k_mean = mean_k,
      gamma = gamma,
      n_vac = n_vac_allocs, n_rep = n_rep,
      timepoints = timepoints, init_I = init_I_nw,
      mc.cores = inner_cores),
    stop("Unknown model: ", job$model))
  setDT(ve)
  ve[, model        := job$model]
  ve[, variant      := job$variant]
  ve[, alpha        := job$alpha]
  ve[, network_seed := job$network_seed %||% NA_integer_]
  ve[]
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

message(glue("Dispatching {length(jobs)} jobs across outer_cores={outer_cores}, inner_cores={inner_cores}"))

results <- parallel::mclapply(jobs, run_job,
                              mc.cores       = outer_cores,
                              mc.preschedule = FALSE)

errs <- vapply(results, inherits, logical(1), what = "try-error")
if (any(errs)) {
  warning(glue("{sum(errs)}/{length(jobs)} jobs failed; dropping."))
  results <- results[!errs]
}
all_res <- rbindlist(results, fill = TRUE)
fwrite(all_res, file.path(out_dir, "raw.csv"))
message(glue("Wrote {out_dir}/raw.csv ({nrow(all_res)} rows)"))

# ---------------------------------------------------------------------------
# Aggregate + plot
# ---------------------------------------------------------------------------

ve_dt <- all_res[method == "full_stoch"]
ve_dt[, VE := 1 - eate]

# Combined label for colour aesthetic (empty variant => just the model)
ve_dt[, model_variant := ifelse(variant == "" | is.na(variant),
                                as.character(model),
                                paste0(model, "_", variant))]

# Median + 95% quantile across all (network_seed, sim) at fixed
# (t, model_variant, alpha). For network this pools over the
# n_network_seeds contact matrices AND the n_vac_allocs allocations
# within each. For frailty it pools over allocations.
summary_dt <- ve_dt[, .(VE_med = median(VE, na.rm = TRUE),
                        VE_lo  = quantile(VE, 0.025, na.rm = TRUE),
                        VE_hi  = quantile(VE, 0.975, na.rm = TRUE),
                        n      = .N),
                    by = .(t, model, variant, model_variant, alpha)]
fwrite(summary_dt, file.path(out_dir, "summary.csv"))

# Plot 1: VE(t), facet by alpha, colour by (model + variant)
p_by_alpha <- ggplot(summary_dt,
                     aes(x = t, y = VE_med,
                         colour = model_variant, fill = model_variant,
                         group  = model_variant)) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              alpha = 0.18, colour = NA) +
  geom_line(size = 1) +
  facet_wrap(~ alpha, labeller = label_both) +
  scale_colour_viridis_d(name = "model") +
  scale_fill_viridis_d(name   = "model") +
  theme_minimal(base_size = 13) +
  labs(x = "t", y = "VE = 1 - EATE (full_stoch)",
       title = glue("VE(t) by model — beta = {beta}, gamma = {gamma}, N = {N}"))
ggsave(file.path(out_dir, "ve_by_alpha.png"),
       p_by_alpha,
       width = 13, height = 8, dpi = 130)

# Plot 2: VE(t), facet by (model + variant), colour by alpha
p_by_model <- ggplot(summary_dt,
                     aes(x = t, y = VE_med,
                         colour = factor(alpha), fill = factor(alpha),
                         group  = interaction(model_variant, alpha))) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              alpha = 0.15, colour = NA) +
  geom_line(size = 1) +
  facet_wrap(~ model_variant, scales = "free_y") +
  scale_colour_viridis_d(name = "alpha") +
  scale_fill_viridis_d(name   = "alpha") +
  theme_minimal(base_size = 13) +
  labs(x = "t", y = "VE = 1 - EATE",
       title = glue("VE(t) per model — beta = {beta}, gamma = {gamma}, N = {N}"))
ggsave(file.path(out_dir, "ve_by_model.png"),
       p_by_model,
       width = 13, height = 9, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
