# VE(t) by vaccine susceptibility alpha, across model types (paper fig 3).
# One facet per alpha in `alphas`; curves:
#   - Linear                 exposure-only, no heterogeneity (stochastic EATE)
#   - Linear, CV^2 = 2       exposure-only with a mean-1 gamma frailty in risk
#                            (shape 1/CV^2), closed form. No transmission, so
#                            no indirect effects: EATE = CIR.
#   - SIR, random mixing     one homogeneously mixing population
#   - Network                Pareto contact graph, exponent `pl_alpha`, one
#                            line per graph (n_network_seeds graphs)
#
# The stochastic EATE functions average over vaccine allocations and n_rep
# dust replicates per allocation.
#
# Extracted from network_alpha_sweep.R. 14 jobs; the network ones take ~3 min
# each, so about 6 min on 8 cores. Respects SLURM_CPUS_PER_TASK.
#
#   Rscript ve_by_alpha.R

library(dplyr)
library(data.table)
library(ggplot2)
library(glue)
library(RColorBrewer)

source("det_model.R")
source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Parameters (edit here)
# ---------------------------------------------------------------------------

alphas       <- c(0.2, 0.5)   # one facet each
pl_alpha     <- 1.5           # Pareto exponent of the network
cv2_linear   <- 2             # frailty CV^2 for the heterogeneous linear model

beta         <- 1.5   # transmission rate for the SIR and network models
                      # (R0 = beta/gamma = 7.5 in the homogeneous limit)
beta_linear  <- 0.1   # linear model: cum FOI = beta_linear * t = 5 at t = 50,
                      # unvac attack rate ~ 1 - exp(-5) ~ 0.99
gamma        <- 0.2
N            <- 500
t_horizon    <- 50
mean_k       <- 6
init_I_2g    <- c(10, 10)   # initial infected for the SIR models
init_I_nw    <- 20

n_network_seeds <- 5     # independent Pareto graphs
n_vac_allocs    <- 10    # vaccine allocations per network EATE
n_rep           <- 1000  # dust replicates per allocation

outer_cores <- {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (nzchar(slurm_cpus)) as.integer(slurm_cpus)
  else max(1L, parallel::detectCores() - 1L)
}
inner_cores <- 1

out_dir <- "output/ve_by_alpha"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

timepoints <- seq(1, t_horizon, 1)

# ---------------------------------------------------------------------------
# Networks and job list
# ---------------------------------------------------------------------------

networks <- lapply(seq_len(n_network_seeds), function(ns) {
  set.seed(1000L * round(pl_alpha * 10) + ns)
  get_conact_matrix_pl(N, alpha = pl_alpha, mean_k = mean_k)
})
set.seed(NULL)

jobs <- list()
for (a in alphas) {
  jobs[[length(jobs) + 1L]] <- list(model = "linear",       alpha = a)
  jobs[[length(jobs) + 1L]] <- list(model = "sir",          alpha = a)
  for (ns in seq_len(n_network_seeds))
    jobs[[length(jobs) + 1L]] <- list(model = "network", alpha = a, network_seed = ns)
}
message(glue("Built {length(jobs)} jobs."))

`%||%` <- function(a, b) if (is.null(a)) b else a

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
    network = get_stoch_eate_network(
      beta = beta, susceptibility = susc, f = 0.5, N = N,
      t = t_horizon, c_ij = networks[[job$network_seed]], k_mean = mean_k,
      gamma = gamma,
      n_vac = n_vac_allocs, n_rep = n_rep,
      timepoints = timepoints, init_I = init_I_nw,
      mc.cores = inner_cores),
    stop("Unknown model: ", job$model))
  setDT(ve)
  ve[, `:=`(model = job$model, alpha = job$alpha,
            network_seed = job$network_seed %||% NA_integer_)]
  ve[]
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

message(glue("Dispatching {length(jobs)} jobs across outer_cores={outer_cores}"))
results <- parallel::mclapply(jobs, run_job, mc.cores = outer_cores,
                              mc.preschedule = FALSE)

errs <- vapply(results, inherits, logical(1), what = "try-error")
if (any(errs)) {
  message(glue("=== {sum(errs)}/{length(jobs)} jobs failed ==="))
  for (idx in which(errs)) {
    job <- jobs[[idx]]
    cond_msg <- attr(results[[idx]], "condition")$message %||% as.character(results[[idx]])
    message(glue("  [{idx}] model={job$model} alpha={job$alpha} ",
                 "net_seed={job$network_seed %||% NA} -> {cond_msg}"))
  }
  results <- results[!errs]
}
all_res <- rbindlist(results, fill = TRUE)
fwrite(all_res, file.path(out_dir, "raw.csv"))

# ---------------------------------------------------------------------------
# Aggregate, add the closed-form heterogeneous linear model, plot
# ---------------------------------------------------------------------------

ve_dt <- all_res[method == "full_stoch"]
ve_dt[, VE := 1 - eate]
ve_dt[, line_id := ifelse(model == "network", sprintf("network_n%02d", network_seed), model)]
summary_dt <- ve_dt[, .(VE = mean(VE, na.rm = TRUE), ar_c = mean(ar_c, na.rm = TRUE),
                        n = .N),
                    by = .(t, model, line_id, alpha)]

# Mean-1 gamma frailty with shape k = 1/CV^2: S(t) = (1 + lambda t / k)^(-k),
# vaccinated at rate alpha * lambda. VE = 1 - (1 - S_v) / (1 - S_c).
k_lin <- 1 / cv2_linear
lin_het <- CJ(t = timepoints, alpha = alphas)
lin_het[, VE := 1 - (1 - (1 + alpha * beta_linear * t / k_lin)^(-k_lin)) /
                    (1 - (1 + beta_linear * t / k_lin)^(-k_lin))]
lin_het[, ar_c := 1 - (1 + beta_linear * t / k_lin)^(-k_lin)]
lin_het[, `:=`(model = "linear_het", line_id = "linear_het", n = NA_integer_)]
summary_dt <- rbind(summary_dt, lin_het, use.names = TRUE)
fwrite(summary_dt, file.path(out_dir, "summary.csv"))

model_labels <- c(
  linear       = "Linear",
  linear_het   = sprintf("Linear, risk CV² = %g", cv2_linear),
  sir          = "SIR, random mixing",
  network      = sprintf("Network (Pareto exp. = %g)", pl_alpha))
summary_dt[, display := factor(model_labels[model], levels = model_labels)]
summary_dt[, alpha_f := factor(alpha, levels = alphas, labels = paste("\u03b1 =", alphas))]

# A: VE against time. B: VE against the control-arm attack rate, which puts
# models with different epidemic speeds on a common depletion scale.
ve_panel <- function(x, xlab) {
  ggplot(summary_dt, aes(x = .data[[x]], y = VE, colour = display, group = line_id)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ alpha_f, scales = "free_y") +
    scale_colour_brewer(NULL, palette = "Dark2") +
    theme_minimal() +
    theme(text = element_text(size = 20), legend.position = "bottom",
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = xlab, y = "VE")
}
p_t  <- ve_panel("t", "t") + theme(legend.position = "none")
p_ar <- ve_panel("ar_c", "Cumulative incidence control") + coord_cartesian(xlim = c(0, 1))
p_by_alpha <- cowplot::plot_grid(p_t, p_ar, ncol = 1, labels = c("A", "B"),
                                 rel_heights = c(1, 1.25))
ggsave("output/ve_by_alpha.png", p_by_alpha, bg = "white",
       width = 13, height = 13, dpi = 130)

message(glue("Done. output/ve_by_alpha.png, data in {out_dir}/"))
