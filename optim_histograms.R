# Stochastic simulator distributions across five model variants:
#   1. linear            — no contact dynamics, run_stoch_linear_dust
#   2. sir               — homogeneous SIR via run_stoch_cd_dust
#   3. sir_sus_frailty   — SIR + susceptibility frailty (sd > 0)
#   4. sir_trans_frailty — SIR + transmissibility frailty (sd_trans > 0)
#   5. network           — Pareto contact network via run_stoch_network
#
# For each model we Nelder-Mead optimise (beta, alpha) so that the
# stochastic means E[C1], E[C2] at t* match the targets (data_C1, data_C2):
#   loss = (mean(C1) - data_C1)^2 + (mean(C2) - data_C2)^2
# Then we run n_sim_final stochastic replicates with the optimised
# parameters and plot histograms of C1, C2, and C1/C2 across replicates.
#
# Usage:
#   Rscript optim_histograms.R

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

N_cont    <- 100
N_vac     <- 100
data_C1   <- 50        # target mean unvac cases at t*
data_C2   <- 25        # target mean vac cases at t*
t_star    <- 8

n_sim_opt   <- 200     # stochastic replicates per loss evaluation
n_sim_final <- 1000    # stochastic replicates for the histogram run

gamma     <- 1         # project convention
I_ini_2g  <- c(10, 10) # seeds per group for SIR / frailty
init_I_nw <- 2         # initial infected count for network
dt        <- 0.01
cores     <- 4

# Frailty params
sd_sus     <- 0.3
sd_trans_v <- 0.3
n_frailty  <- 10

# Network params
pl_alpha     <- 3
mean_k       <- 6
network_seed <- 1
alloc_seed   <- 1

# Optimiser starting point and iterations
beta_ini    <- 1.0
alpha_ini   <- 0.5
optim_maxit <- 80

out_dir <- "output/optim_histograms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Build the network and vac allocation once so optimisation/finals share them
# ---------------------------------------------------------------------------

N_total <- N_cont + N_vac

set.seed(network_seed)
c_ij_fixed <- get_conact_matrix_pl(N_total, alpha = pl_alpha, mean_k = mean_k)
set.seed(NULL)

set.seed(alloc_seed)
vac_fixed <- sample(seq_len(N_total), N_vac)
set.seed(NULL)

# ---------------------------------------------------------------------------
# Model wrappers: each returns a data.table at time == t_star with C1, C2
# ---------------------------------------------------------------------------

at_tstar <- function(out) {
  setDT(out)
  out[time == t_star, .(sim, C1, C2)]
}

run_linear <- function(beta, alpha, n_sim) {
  out <- run_stoch_linear_dust(beta = beta,
                               N = c(N_cont, N_vac),
                               susceptibility = c(1, alpha),
                               t = t_star, dt = dt,
                               timepoints = seq(1, t_star, 1),
                               n_sim = n_sim, cores = cores)
  # linear output has columns S1, S2, C1, C2 — rename C1/C2 stays as is
  at_tstar(out)
}

run_sir <- function(beta, alpha, n_sim) {
  out <- run_stoch_cd_dust(matrix(rep(1, 4), nrow = 2),
                           beta = beta, N = c(N_cont, N_vac),
                           t = t_star, I_ini = I_ini_2g,
                           susceptibility = c(1, alpha),
                           gamma = gamma, dt = dt,
                           timepoints = seq(1, t_star, 1),
                           n_sim = n_sim, cores = cores)
  at_tstar(out)
}

run_sir_sus_frailty <- function(beta, alpha, n_sim) {
  out <- run_stoch_frailty_cd(sd = sd_sus, sd_trans = 0,
                              susceptibility = c(1, alpha),
                              beta = beta,
                              N = N_total / 2,                # 2*N convention
                              t = t_star, n_frailty = n_frailty,
                              gamma = gamma, I_ini_total = sum(I_ini_2g),
                              timepoints = seq(1, t_star, 1),
                              n_sim = n_sim, cores = cores,
                              method = "dust", dt = dt, f = 0.5)
  at_tstar(out)
}

run_sir_trans_frailty <- function(beta, alpha, n_sim) {
  out <- run_stoch_frailty_cd(sd = 0, sd_trans = sd_trans_v,
                              susceptibility = c(1, alpha),
                              beta = beta,
                              N = N_total / 2,
                              t = t_star, n_frailty = n_frailty,
                              gamma = gamma, I_ini_total = sum(I_ini_2g),
                              timepoints = seq(1, t_star, 1),
                              n_sim = n_sim, cores = cores,
                              method = "dust", dt = dt, f = 0.5)
  at_tstar(out)
}

run_network <- function(beta, alpha, n_sim) {
  out <- run_stoch_network(beta = beta, N = N_total,
                           susceptibility = c(1, alpha),
                           t = t_star, c_ij = c_ij_fixed, vac = vac_fixed,
                           k_mean = mean_k, gamma = gamma,
                           dt = dt, timepoints = seq(1, t_star, 1),
                           n_sim = n_sim, cores = cores, I_ini = init_I_nw)
  at_tstar(out)
}

models <- list(
  linear            = run_linear,
  sir               = run_sir,
  sir_sus_frailty   = run_sir_sus_frailty,
  sir_trans_frailty = run_sir_trans_frailty,
  network           = run_network
)

# ---------------------------------------------------------------------------
# Optimise (beta, alpha) per model on log scale so positivity is automatic
# ---------------------------------------------------------------------------

make_loss <- function(model_fn) {
  function(log_par) {
    par   <- exp(log_par)
    beta  <- par[1]
    alpha <- par[2]
    out <- tryCatch(model_fn(beta, alpha, n_sim_opt),
                    error = function(e) NULL)
    if (is.null(out) || nrow(out) == 0) return(1e9)
    (mean(out$C1) - data_C1)^2 + (mean(out$C2) - data_C2)^2
  }
}

opt_results <- list()
for (mname in names(models)) {
  message(glue("Optimising {mname}..."))
  loss <- make_loss(models[[mname]])
  o <- optim(par = log(c(beta_ini, alpha_ini)), fn = loss,
             method = "Nelder-Mead",
             control = list(maxit = optim_maxit, reltol = 1e-3))
  pars <- exp(o$par)
  opt_results[[mname]] <- list(beta = pars[1], alpha = pars[2],
                               loss = o$value, convergence = o$convergence)
  message(glue("  beta = {round(pars[1], 4)}  alpha = {round(pars[2], 4)}  ",
               "loss = {round(o$value, 3)}  conv = {o$convergence}"))
}

opt_dt <- rbindlist(lapply(names(opt_results), function(m) {
  r <- opt_results[[m]]
  data.table(model = m, beta = r$beta, alpha = r$alpha,
             loss = r$loss, convergence = r$convergence)
}))
fwrite(opt_dt, file.path(out_dir, "optimised_params.csv"))
message("\nOptimised parameters:")
print(opt_dt)

# ---------------------------------------------------------------------------
# Final stochastic simulations + histograms
# ---------------------------------------------------------------------------

final_results <- list()
for (mname in names(models)) {
  r <- opt_results[[mname]]
  message(glue("\nFinal simulation for {mname} (beta={round(r$beta, 4)}, alpha={round(r$alpha, 4)})..."))
  out <- models[[mname]](r$beta, r$alpha, n_sim_final)
  out[, model := mname]
  final_results[[mname]] <- out
}

all_res <- rbindlist(final_results)
all_res[, ratio := C1 / C2]
all_res[!is.finite(ratio), ratio := NA_real_]

fwrite(all_res, file.path(out_dir, "final_simulations.csv"))

# Preserve the model order on the facets
all_res[, model := factor(model, levels = names(models))]
target_ratio <- data_C1 / data_C2

plot_hist <- function(df, var, title, xlab, target, fill) {
  ggplot(df, aes(x = .data[[var]])) +
    geom_histogram(bins = 40, fill = fill, alpha = 0.75) +
    geom_vline(xintercept = target, colour = "red", linetype = "dashed") +
    facet_wrap(~ model, scales = "free_y", ncol = 1) +
    theme_minimal(base_size = 13) +
    labs(x = xlab, y = "count", title = title)
}

ggsave(file.path(out_dir, "hist_C1.png"),
       plot_hist(all_res, "C1",
                 glue("C1 distribution at t={t_star} (target = {data_C1})"),
                 "C1 (unvac cases)", data_C1, "steelblue"),
       width = 8, height = 12, dpi = 130)

ggsave(file.path(out_dir, "hist_C2.png"),
       plot_hist(all_res, "C2",
                 glue("C2 distribution at t={t_star} (target = {data_C2})"),
                 "C2 (vac cases)", data_C2, "darkorange"),
       width = 8, height = 12, dpi = 130)

ggsave(file.path(out_dir, "hist_ratio.png"),
       plot_hist(all_res[!is.na(ratio)], "ratio",
                 glue("C1 / C2 at t={t_star} (target = {round(target_ratio, 3)})"),
                 "C1 / C2", target_ratio, "purple"),
       width = 8, height = 12, dpi = 130)

message(glue("\nDone. Outputs in {out_dir}/"))
