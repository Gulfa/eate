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
beta         <- 1.5   # transmission rate for sir / frailty / network
                      # (R0 = beta/gamma in the homogeneous limit).
beta_linear  <- 0.1   # linear model has no recovery / herd effect, so
                      # cum FOI = beta * t. beta_linear = 0.1 over t = 8
                      # gives an unvac attack rate ~ 1 - exp(-0.8) ~ 0.55,
                      # comparable to R0 = 1.5 SIR's final size.
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
      beta = beta_linear, susceptibility = susc, f = 0.5, N = N,
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

# Two labels:
#   color_group: what determines line colour. For network this is
#     model + variant WITHOUT network_seed, so all draws of the same
#     pl_alpha share a colour.
#   line_id: uniquely identifies a curve. For network this includes
#     network_seed so each draw is its own line.
ve_dt[, color_group := ifelse(
  variant == "" | is.na(variant),
  as.character(model),
  paste0(model, "_", variant))]
ve_dt[, line_id := ifelse(
  model == "network",
  sprintf("%s_n%02d", color_group, network_seed),
  color_group)]

# Mean VE at fixed (t, line_id, alpha). Network: mean over inner
# n_vac_allocs allocations only, different graphs stay separate.
# Frailty: mean over allocations.
summary_dt <- ve_dt[, .(VE = mean(VE, na.rm = TRUE), n = .N),
                    by = .(t, model, variant, color_group, line_id, alpha)]
fwrite(summary_dt, file.path(out_dir, "summary.csv"))

# Plot 1: VE(t), facet by alpha. Colour by color_group so network
# draws of the same pl_alpha share a colour; group by line_id so each
# network seed is its own line.
p_by_alpha <- ggplot(summary_dt,
                     aes(x = t, y = VE,
                         colour = color_group,
                         group  = line_id)) +
  geom_line(size = 1) +
  facet_wrap(~ alpha, labeller = label_both) +
  scale_colour_viridis_d(name = "model") +
  theme_minimal(base_size = 13) +
  labs(x = "t", y = "VE = 1 - EATE (mean over allocations)",
       title = glue("VE(t) by model — beta = {beta} (linear beta = {beta_linear}), gamma = {gamma}, N = {N}"))
ggsave(file.path(out_dir, "ve_by_alpha.png"),
       p_by_alpha,
       width = 13, height = 8, dpi = 130)

# Plot 2: VE(t), facet by color_group, colour by alpha. Multiple
# network draws of the same pl_alpha share the same facet.
p_by_model <- ggplot(summary_dt,
                     aes(x = t, y = VE,
                         colour = factor(alpha),
                         group  = interaction(line_id, alpha))) +
  geom_line(size = 1) +
  facet_wrap(~ color_group, scales = "free_y") +
  scale_colour_viridis_d(name = "alpha") +
  theme_minimal(base_size = 13) +
  labs(x = "t", y = "VE (mean over allocations)",
       title = glue("VE(t) per model — beta = {beta} (linear beta = {beta_linear}), gamma = {gamma}, N = {N}"))
ggsave(file.path(out_dir, "ve_by_model.png"),
       p_by_model,
       width = 13, height = 9, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
