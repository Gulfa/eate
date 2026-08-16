# VE + AVE with uncertainty from a supplied (beta_hat, alpha_hat) posterior
# for a homogeneous SIR model, compared against the naive difference-of-
# binomials interval a trial would compute from the observed C1, C2.
#
# The point: propagating parameter uncertainty through the SIR (posterior
# sample -> counterfactual EATE + AVE per sample -> quantile) gives a
# principled AVE distribution that accounts for both stochastic case
# generation AND the parameter uncertainty implied by the fit. The Wald
# difference-of-binomials interval, by contrast, only accounts for
# residual sampling noise in the observed C1, C2 (conditional on a fixed
# underlying AR), and typically underestimates the true uncertainty when
# the fit is loose.
#
# Usage: Rscript ave_from_sir_fit.R  (edit inputs at the top)

library(dplyr)
library(data.table)
library(ggplot2)
library(glue)
library(RColorBrewer)

source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Inputs — edit here
# ---------------------------------------------------------------------------

# Fitted point estimate
beta_hat   <- 1.5
alpha_hat  <- 0.5

# Posterior covariance in (beta, alpha) order. Two ways to specify:
#   - Explicit 2x2 matrix, or
#   - marginal SDs + correlation.
sd_beta    <- 0.10
sd_alpha   <- 0.06
cor_ba     <- -0.2
cov_hat    <- matrix(c(sd_beta^2,             cor_ba * sd_beta * sd_alpha,
                       cor_ba * sd_beta * sd_alpha, sd_alpha^2),
                     nrow = 2L, byrow = TRUE)

# Study / data design that produced this fit
N_cont     <- 100
N_vac      <- 100
data_C1    <- 50     # observed control cases at t_star
data_C2    <- 25     # observed vac      cases at t_star
t_star     <- 8
gamma      <- 1
I_ini_2g   <- c(10, 10)
dt         <- 0.01

# EATE / AVE knobs
K          <- 200    # posterior parameter samples
n_vac      <- 5      # inner vac allocations per EATE call
n_rep      <- 200    # dust replicates per allocation
cores      <- 4
seed_grid  <- 1234L

out_dir <- "output/ave_from_sir_fit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Posterior samples
# ---------------------------------------------------------------------------

set.seed(seed_grid)
L <- tryCatch(chol(cov_hat),
              error = function(e) chol(cov_hat + diag(1e-10, 2L)))
Z <- matrix(rnorm(2L * K), nrow = K, ncol = 2L)
draws <- sweep(Z %*% L, 2L, c(beta_hat, alpha_hat), "+")
colnames(draws) <- c("beta", "alpha")
# Keep only samples inside the positive orthant.
keep  <- draws[, "beta"] > 0 & draws[, "alpha"] > 0
draws <- draws[keep, , drop = FALSE]
set.seed(NULL)
message(glue("Using {nrow(draws)} / {K} posterior samples (positive orthant)."))

# ---------------------------------------------------------------------------
# Run get_stoch_eate_sir at each sample
# ---------------------------------------------------------------------------

N_total <- N_cont + N_vac
vac_frac <- N_vac / N_total
timepoints <- seq_len(t_star)

rows <- vector("list", nrow(draws))
for (k in seq_len(nrow(draws))) {
  bk <- draws[k, "beta"]; ak <- draws[k, "alpha"]
  if (k %% 25L == 0L || k == 1L)
    message(glue("  posterior sample {k}/{nrow(draws)}"))
  ve <- get_stoch_eate_sir(
    beta = bk, susceptibility = c(1, ak), f = vac_frac, N = N_total,
    t = t_star, gamma = gamma, I_ini = I_ini_2g,
    n_vac = n_vac, n_rep = n_rep,
    dt = dt, timepoints = timepoints,
    mc.cores = cores)
  setDT(ve)
  ve[, param_sample := k]
  rows[[k]] <- ve
}
ve_all <- rbindlist(rows)

# ---------------------------------------------------------------------------
# Aggregate at t = t_star, method = "full_stoch"
# ---------------------------------------------------------------------------

final <- ve_all[t == t_star & method == "full_stoch"]
final[, VE  := 1 - eate]
# AVE column is populated by get_stoch_eate_sir as (denom - num) / N_total,
# so positive means the vaccine averted cases.
setnames(final, "ave", "AVE")

# Quantiles over param_sample x sim (posterior + Monte-Carlo pooled)
ve_summary <- quantile(final$VE,  probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
ave_summary <- quantile(final$AVE, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# ---------------------------------------------------------------------------
# Naive difference-of-binomials interval (Wald)
# ---------------------------------------------------------------------------

ar_c_obs <- data_C1 / N_cont
ar_v_obs <- data_C2 / N_vac
ave_naive <- ar_c_obs - ar_v_obs
se_naive  <- sqrt(ar_c_obs * (1 - ar_c_obs) / N_cont +
                  ar_v_obs * (1 - ar_v_obs) / N_vac)
ave_naive_lo <- ave_naive - 1.96 * se_naive
ave_naive_hi <- ave_naive + 1.96 * se_naive

# Corresponding naive VE (1 - AR_v/AR_c) with delta-method SE, for reference.
ve_naive     <- 1 - ar_v_obs / ar_c_obs
# delta method on log-risk-ratio scale:
#   Var(log(RR)) ≈ (1 - p_v)/(n_v * p_v) + (1 - p_c)/(n_c * p_c)
log_rr_se <- sqrt((1 - ar_v_obs) / (N_vac * ar_v_obs) +
                  (1 - ar_c_obs) / (N_cont * ar_c_obs))
rr <- ar_v_obs / ar_c_obs
rr_lo <- rr * exp(-1.96 * log_rr_se)
rr_hi <- rr * exp( 1.96 * log_rr_se)
ve_naive_lo <- 1 - rr_hi
ve_naive_hi <- 1 - rr_lo

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

message("\n===================================================================")
message(glue("Fitted point: beta = {beta_hat}, alpha = {alpha_hat}"))
message(glue("Posterior sd: sd_beta = {sd_beta}, sd_alpha = {sd_alpha}, ",
             "cor = {cor_ba}"))
message(glue("Study: N_cont = {N_cont}, N_vac = {N_vac}, ",
             "C1 = {data_C1}, C2 = {data_C2}, t* = {t_star}"))
message("--- Posterior-propagated (via SIR + K samples) ---")
message(glue("  VE  at t*: median = {round(ve_summary[2], 3)}, ",
             "95%% CI = [{round(ve_summary[1], 3)}, {round(ve_summary[3], 3)}]"))
message(glue("  AVE at t*: median = {round(ave_summary[2], 3)}, ",
             "95%% CI = [{round(ave_summary[1], 3)}, {round(ave_summary[3], 3)}]"))
message("--- Naive difference-of-binomials (Wald) ---")
message(glue("  VE  (log-RR delta): {round(ve_naive, 3)}, ",
             "95%% CI = [{round(ve_naive_lo, 3)}, {round(ve_naive_hi, 3)}]"))
message(glue("  AVE (Wald)        : {round(ave_naive, 3)}, ",
             "95%% CI = [{round(ave_naive_lo, 3)}, {round(ave_naive_hi, 3)}]"))
message("===================================================================\n")

fwrite(final, file.path(out_dir, "posterior_ve_ave_final.csv"))

# ---------------------------------------------------------------------------
# Plot: AVE posterior histogram with both CIs overlaid
# ---------------------------------------------------------------------------

p_ave <- ggplot(final, aes(x = AVE)) +
  geom_histogram(bins = 40, fill = brewer.pal(3L, "Dark2")[1], alpha = 0.7) +
  # Posterior CI (blue)
  geom_vline(xintercept = ave_summary[2], colour = brewer.pal(3L, "Dark2")[1], size = 1) +
  geom_vline(xintercept = c(ave_summary[1], ave_summary[3]),
             colour = brewer.pal(3L, "Dark2")[1], linetype = "dashed") +
  # Naive Wald CI (orange)
  geom_vline(xintercept = ave_naive,   colour = brewer.pal(3L, "Dark2")[2], size = 1) +
  geom_vline(xintercept = c(ave_naive_lo, ave_naive_hi),
             colour = brewer.pal(3L, "Dark2")[2], linetype = "dashed") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "AVE (per-person absolute effect at t*)",
       y = "count",
       title = glue("AVE at t = {t_star}: posterior propagation (blue) vs Wald difference of binomials (orange)"))
ggsave(file.path(out_dir, "ave_posterior_vs_wald.png"), p_ave,
       width = 10, height = 5, dpi = 130)

# Same for VE, for completeness
p_ve <- ggplot(final, aes(x = VE)) +
  geom_histogram(bins = 40, fill = brewer.pal(3L, "Dark2")[1], alpha = 0.7) +
  geom_vline(xintercept = ve_summary[2], colour = brewer.pal(3L, "Dark2")[1], size = 1) +
  geom_vline(xintercept = c(ve_summary[1], ve_summary[3]),
             colour = brewer.pal(3L, "Dark2")[1], linetype = "dashed") +
  geom_vline(xintercept = ve_naive,      colour = brewer.pal(3L, "Dark2")[2], size = 1) +
  geom_vline(xintercept = c(ve_naive_lo, ve_naive_hi),
             colour = brewer.pal(3L, "Dark2")[2], linetype = "dashed") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "VE at t*", y = "count",
       title = glue("VE at t = {t_star}: posterior propagation (blue) vs Wald log-RR (orange)"))
ggsave(file.path(out_dir, "ve_posterior_vs_wald.png"), p_ve,
       width = 10, height = 5, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
