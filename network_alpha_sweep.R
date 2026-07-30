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
library(RColorBrewer)

source("det_model.R")
source("stoch_model.R")
source("utils.R")

# Dark2 palette extended to n groups via colorRampPalette when needed.
dark2_scale <- function(n, name = NULL) {
  cols <- if (n <= 8L) brewer.pal(max(n, 3L), "Dark2")[seq_len(n)]
          else colorRampPalette(brewer.pal(8L, "Dark2"))(n)
  scale_colour_manual(name = name, values = cols)
}

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
  message(glue("=== {sum(errs)}/{length(jobs)} jobs failed ==="))
  for (idx in which(errs)) {
    job <- jobs[[idx]]
    cond_msg <- attr(results[[idx]], "condition")$message
    if (is.null(cond_msg)) cond_msg <- as.character(results[[idx]])
    message(glue("  [{idx}] model={job$model}  variant={job$variant}  ",
                 "alpha={job$alpha}  net_seed={job$network_seed %||% NA}  ",
                 "-> {cond_msg}"))
  }
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

n_groups   <- length(unique(summary_dt$color_group))
n_alphas   <- length(unique(summary_dt$alpha))

# Publication-ready display name mapping. Ordered so the legend flows
# from simplest (linear) to most structured (network).
display_name <- function(color_group) {
  out <- character(length(color_group))
  for (i in seq_along(color_group)) {
    g <- color_group[i]
    out[i] <- if (g == "linear") "Linear"
    else if (g == "sir") "SIR (homogeneous)"
    else if (grepl("^sir_sus_frailty_sd", g)) {
      sd <- sub("^sir_sus_frailty_sd", "", g)
      sprintf("SIR + sus. frailty (σ = %s)", sd)
    } else if (grepl("^sir_trans_frailty_sd", g)) {
      sd <- sub("^sir_trans_frailty_sd", "", g)
      sprintf("SIR + trans. frailty (σ = %s)", sd)
    } else if (grepl("^network_pa", g)) {
      pa <- sub("^network_pa", "", g)
      sprintf("Network (Pareto exp. = %s)", pa)
    } else g
  }
  out
}

# Build an ordered factor so legend + colour cycle in a sensible order.
groups_ordered <- unique(summary_dt$color_group)
order_score <- function(g) {
  if (g == "linear")                          return(1)
  if (g == "sir")                             return(2)
  if (grepl("^sir_sus_frailty_sd", g))        return(3 + as.numeric(sub("^sir_sus_frailty_sd", "", g)))
  if (grepl("^sir_trans_frailty_sd", g))      return(5 + as.numeric(sub("^sir_trans_frailty_sd", "", g)))
  if (grepl("^network_pa", g))                return(7 + as.numeric(sub("^network_pa", "", g)))
  99
}
groups_ordered <- groups_ordered[order(sapply(groups_ordered, order_score))]
group_labels   <- setNames(display_name(groups_ordered), groups_ordered)

summary_dt[, display := factor(color_group, levels = groups_ordered,
                                labels = group_labels[groups_ordered])]

# Plot 1: VE(t), facet by alpha. Publication style — Dark2 palette,
# white background, legend at bottom, larger text.
p_by_alpha <- ggplot(summary_dt,
                     aes(x = t, y = VE,
                         colour = display,
                         group  = line_id)) +
  geom_line(size = 0.9) +
  facet_wrap(~ alpha, labeller = label_both, scales = "free_y") +
  dark2_scale(n_groups, name = NULL) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "t", y = "VE")
ggsave(file.path(out_dir, "ve_by_alpha.png"),
       p_by_alpha,
       width = 13, height = 8, dpi = 130)

# Plot 2: VE(t), facet by display (nice model label), colour by alpha.
# Multiple network draws of the same pl_alpha share the same facet.
p_by_model <- ggplot(summary_dt,
                     aes(x = t, y = VE,
                         colour = factor(alpha),
                         group  = interaction(line_id, alpha))) +
  geom_line(size = 0.9) +
  facet_wrap(~ display, scales = "free_y") +
  dark2_scale(n_alphas, name = "alpha") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "t", y = "VE")
ggsave(file.path(out_dir, "ve_by_model.png"),
       p_by_model,
       width = 13, height = 9, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
