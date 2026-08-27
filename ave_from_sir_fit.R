# VE + AVE with uncertainty from a supplied (beta_hat, alpha_hat) posterior,
# compared against the naive difference-of-binomials interval a trial would
# compute from the observed C1, C2. Runs the diagnostic for two model
# families:
#
#   - "linear" : exposure-only, no interaction. Closed-form MLE from
#                (C1, C2), covariance via the sandwich on the linear
#                simulator. Should match Wald exactly (2 params, 2
#                moments, bijection between (beta, alpha) and (C1, C2)).
#   - "sir"    : homogeneous SIR. User supplies the fit and its
#                covariance (from a real fit / run_fit_array). Should
#                produce wider CIs than Wald because Var(C_i | theta) is
#                super-binomial and J has depletion-driven shrinkage.
#
# The intended sanity check: linear = Wald ✓, SIR ⩾ Wald ✓. If linear
# also comes out wider than Wald, something is still off upstream
# (posterior sampling, sandwich derivation, or the AVE column).
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

# Which models to run
models_to_run <- c("linear", "sir")

# --- Shared study / data design ---
N_cont     <- 500
N_vac      <- 500
data_C1    <- 200        # observed control cases at t_star
data_C2    <- 100        # observed vac      cases at t_star
t_star     <- 4
dt         <- 0.01

# --- SIR-specific: supplied fit ---
sir_beta_hat   <- 2.2
sir_alpha_hat  <- 0.49
sir_sd_beta    <- 0.22
sir_sd_alpha   <- 0.06
sir_cor_ba     <- -0.5
sir_gamma      <- 1
sir_I_ini_2g   <- c(10, 10)

# --- Linear-specific: supplied fit (mirrors SIR). Set any of lin_beta_hat /
#     lin_alpha_hat / lin_sd_beta / lin_sd_alpha / lin_cor_ba to NA to fall
#     back to defaults derived from the observed (data_C1, data_C2):
#       - point estimate: closed-form MLE
#           beta  = -log(1 - data_C1/N_cont) / t_star
#           alpha =  log(1 - data_C2/N_vac) / log(1 - data_C1/N_cont)
#       - covariance: sandwich J^-1 Sigma J^-T from the linear simulator
#         (see estimate_posterior_cov). This is the "matches Wald exactly"
#         reference — override only when you want to test a specific fit.
lin_beta_hat   <- 0.13
lin_alpha_hat  <- 0.44
lin_sd_beta    <- 0.009
lin_sd_alpha   <- 0.054
lin_cor_ba     <- -0.6

# --- Diagnostic knobs ---
K              <- 200      # posterior parameter samples
n_vac          <- 5        # inner vac allocations per EATE call
n_rep          <- 200      # dust replicates per allocation
cores          <- 1
seed_grid      <- 1234L
post_cov_n_sim <- 1000     # simulator draws for estimate_posterior_cov (linear)
post_cov_seed  <- 4321L

out_dir <- "output/ave_from_sir_fit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_total  <- N_cont + N_vac        
vac_frac <- N_vac / N_total
timepoints <- seq_len(t_star)

# ---------------------------------------------------------------------------
# Naive difference-of-binomials interval (Wald) — same for both models,
# depends only on the observed data.
# ---------------------------------------------------------------------------

ar_c_obs     <- data_C1 / N_cont
ar_v_obs     <- data_C2 / N_vac
ave_naive    <- ar_c_obs - ar_v_obs
se_naive     <- sqrt(ar_c_obs * (1 - ar_c_obs) / N_cont +
                     ar_v_obs * (1 - ar_v_obs) / N_vac)
ave_naive_lo <- ave_naive - 1.96 * se_naive
ave_naive_hi <- ave_naive + 1.96 * se_naive

ve_naive     <- 1 - ar_v_obs / ar_c_obs
log_rr_se    <- sqrt((1 - ar_v_obs) / (N_vac * ar_v_obs) +
                     (1 - ar_c_obs) / (N_cont * ar_c_obs))
rr           <- ar_v_obs / ar_c_obs
rr_lo        <- rr * exp(-1.96 * log_rr_se)
rr_hi        <- rr * exp( 1.96 * log_rr_se)
ve_naive_lo  <- 1 - rr_hi
ve_naive_hi  <- 1 - rr_lo

# ---------------------------------------------------------------------------
# Linear closed-form MLE from data + sandwich covariance via the simulator
# ---------------------------------------------------------------------------

# For linear:
#   p_c = 1 - exp(-beta * t),  p_v = 1 - exp(-alpha * beta * t)
# so given p_c = data_C1/N_c and p_v = data_C2/N_v:
#   beta  = -log(1 - p_c) / t
#   alpha = log(1 - p_v) / log(1 - p_c)
# Closed-form MLE defaults from the data — used as fallback when
# lin_beta_hat / lin_alpha_hat above are NA.
lin_beta_hat_default  <- -log(1 - ar_c_obs) / t_star
lin_alpha_hat_default <-  log(1 - ar_v_obs) /  log(1 - ar_c_obs)

# Wrap get_stoch_eate_linear behind a (beta, alpha, n_sim, seed) interface
# so estimate_posterior_cov can drive it.
sim_linear <- function(beta, alpha, n_sim, seed = NULL) {
  raw <- run_stoch_linear_dust(
    beta = beta, N = c(N_cont, N_vac),
    susceptibility = c(1, alpha),
    t = t_star, dt = dt, timepoints = timepoints,
    n_sim = n_sim, cores = cores, seed = seed)
  setDT(raw)
  raw[time == t_star, .(sim, C1, C2)]
}

# ---------------------------------------------------------------------------
# Per-model diagnostic runner
# ---------------------------------------------------------------------------

run_diagnostic <- function(model_type) {
  message(glue("\n=== Running diagnostic for model_type = {model_type} ==="))

  if (model_type == "sir") {
    beta_hat  <- sir_beta_hat
    alpha_hat <- sir_alpha_hat
    cov_hat   <- matrix(c(sir_sd_beta^2,
                          sir_cor_ba * sir_sd_beta * sir_sd_alpha,
                          sir_cor_ba * sir_sd_beta * sir_sd_alpha,
                          sir_sd_alpha^2),
                        nrow = 2L, byrow = TRUE)
    sim_fn <- function(beta, alpha, n_sim, seed = NULL) {
      raw <- run_stoch_cd_dust(
        matrix(rep(1, 4), nrow = 2),
        beta = beta, N = c(N_cont, N_vac),
        t = t_star, I_ini = sir_I_ini_2g,
        susceptibility = c(1, alpha),
        gamma = sir_gamma, dt = dt,
        timepoints = timepoints,
        n_sim = n_sim, cores = cores, seed = seed)
      setDT(raw); raw[time == t_star, .(sim, C1, C2)]
    }
    ve_fn <- function(bk, ak) {
      get_stoch_eate_sir(
        beta = bk, susceptibility = c(1, ak), f = vac_frac, N = N_total,
        t = t_star, gamma = sir_gamma, I_ini = sir_I_ini_2g,
        n_vac = n_vac, n_rep = n_rep,
        dt = dt, timepoints = timepoints, mc.cores = cores)
    }
  } else if (model_type == "linear") {
    # Point estimate: use supplied lin_* if given, else closed-form MLE.
    if (is.na(lin_beta_hat) || is.na(lin_alpha_hat)) {
      beta_hat  <- lin_beta_hat_default
      alpha_hat <- lin_alpha_hat_default
      message(glue("  linear point: closed-form MLE from data — ",
                   "beta = {round(beta_hat, 5)}, alpha = {round(alpha_hat, 5)}"))
    } else {
      beta_hat  <- lin_beta_hat
      alpha_hat <- lin_alpha_hat
      message(glue("  linear point: user-supplied — ",
                   "beta = {round(beta_hat, 5)}, alpha = {round(alpha_hat, 5)}"))
    }

    # Covariance: use supplied lin_sd_* / lin_cor_ba if given, else
    # sandwich J^-1 Sigma J^-T from estimate_posterior_cov.
    if (is.na(lin_sd_beta) || is.na(lin_sd_alpha) || is.na(lin_cor_ba)) {
      pcov <- estimate_posterior_cov(sim_linear,
                                     beta = beta_hat, alpha = alpha_hat,
                                     n_sim = post_cov_n_sim,
                                     seed  = post_cov_seed)
      cov_hat <- pcov$cov
      message("  linear cov: sandwich estimate (default)")
    } else {
      cov_hat <- matrix(c(lin_sd_beta^2,
                          lin_cor_ba * lin_sd_beta * lin_sd_alpha,
                          lin_cor_ba * lin_sd_beta * lin_sd_alpha,
                          lin_sd_alpha^2),
                        nrow = 2L, byrow = TRUE)
      message("  linear cov: user-supplied")
    }
    sim_fn  <- sim_linear
    ve_fn <- function(bk, ak) {
      get_stoch_eate_linear(
        beta = bk, susceptibility = c(1, ak), f = vac_frac, N = N_total,
        t = t_star, n_vac = n_vac, n_rep = n_rep,
        dt = dt, timepoints = timepoints, mc.cores = cores)
    }
  } else stop("Unknown model_type: ", model_type)

  message(glue("  fit: beta = {round(beta_hat, 4)}, alpha = {round(alpha_hat, 4)}"))
  message(glue("  sd_beta = {round(sqrt(cov_hat[1,1]), 5)}, ",
               "sd_alpha = {round(sqrt(cov_hat[2,2]), 5)}, ",
               "cor = {round(cov_hat[1,2] / sqrt(cov_hat[1,1] * cov_hat[2,2]), 3)}"))

  # --- Posterior samples ---
  set.seed(seed_grid)
  L <- tryCatch(chol(cov_hat),
                error = function(e) chol(cov_hat + diag(1e-10, 2L)))
  Z <- matrix(rnorm(2L * K), nrow = K, ncol = 2L)
  draws <- sweep(Z %*% L, 2L, c(beta_hat, alpha_hat), "+")
  colnames(draws) <- c("beta", "alpha")
  keep  <- draws[, "beta"] > 0 & draws[, "alpha"] > 0
  draws <- draws[keep, , drop = FALSE]
  set.seed(NULL)
  message(glue("  {nrow(draws)} / {K} MVN samples in positive orthant"))

  # --- Simulate at each posterior sample ---
  rows <- vector("list", nrow(draws))
  for (k in seq_len(nrow(draws))) {
    if (k %% 25L == 0L || k == 1L)
      message(glue("    sample {k}/{nrow(draws)}"))
    ve <- ve_fn(draws[k, "beta"], draws[k, "alpha"])
    setDT(ve); ve[, param_sample := k]
    rows[[k]] <- ve
  }
  ve_all <- rbindlist(rows)

  final <- ve_all[t == t_star & method == "full_stoch"]
  final[, VE := 1 - eate]
  setnames(final, "ave", "AVE")

  # SD-based CI: mean ± 1.96 * SD(pooled). Half the finite-K noise of an
  # empirical 2.5/97.5% quantile at the same K, and consistent with the
  # Laplace / MVN posterior assumption we're using anyway. Named lo/med/hi
  # so downstream ref_post code doesn't need to change.
  ve_ci <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x,   na.rm = TRUE)
    c(m - 1.96 * s, m, m + 1.96 * s)
  }
  ve_q  <- ve_ci(final$VE)
  ave_q <- ve_ci(final$AVE)

  final[, model_type := model_type]
  list(final = final, ve_q = ve_q, ave_q = ave_q,
       beta_hat = beta_hat, alpha_hat = alpha_hat,
       cov_hat = cov_hat)
}

# ---------------------------------------------------------------------------
# Run each model
# ---------------------------------------------------------------------------

diag_results <- setNames(lapply(models_to_run, run_diagnostic), models_to_run)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

message("\n===================================================================")
message(glue("Study: N_cont = {N_cont}, N_vac = {N_vac}, ",
             "C1 = {data_C1}, C2 = {data_C2}, t* = {t_star}"))
message("--- Naive difference-of-binomials (Wald) ---")
message(glue("  VE  (log-RR delta): {round(ve_naive, 3)}, ",
             "95% CI = [{round(ve_naive_lo, 3)}, {round(ve_naive_hi, 3)}] ",
             "half-width = {round((ve_naive_hi - ve_naive_lo) / 2, 3)}"))
message(glue("  AVE (Wald)        : {round(ave_naive, 3)}, ",
             "95% CI = [{round(ave_naive_lo, 3)}, {round(ave_naive_hi, 3)}] ",
             "half-width = {round((ave_naive_hi - ave_naive_lo) / 2, 3)}"))
for (mt in names(diag_results)) {
  r <- diag_results[[mt]]
  message(glue("--- Posterior-propagated: {mt} ---"))
  message(glue("  VE  median = {round(r$ve_q[2], 3)}, ",
               "95% CI = [{round(r$ve_q[1], 3)}, {round(r$ve_q[3], 3)}] ",
               "half-width = {round((r$ve_q[3] - r$ve_q[1]) / 2, 3)}"))
  message(glue("  AVE median = {round(r$ave_q[2], 3)}, ",
               "95% CI = [{round(r$ave_q[1], 3)}, {round(r$ave_q[3], 3)}] ",
               "half-width = {round((r$ave_q[3] - r$ave_q[1]) / 2, 3)}"))
}
message("===================================================================\n")

all_final <- rbindlist(lapply(diag_results, `[[`, "final"))
fwrite(all_final, file.path(out_dir, "posterior_ve_ave_final.csv"))

# ---------------------------------------------------------------------------
# Combined plot: rows = metric (AVE, VE), columns = model_type
# ---------------------------------------------------------------------------

pal_post  <- brewer.pal(3L, "Dark2")[1]
pal_naive <- brewer.pal(3L, "Dark2")[2]

# Data long-format
long_dt <- rbind(
  all_final[, .(model_type, value = AVE, metric = "AVE")],
  all_final[, .(model_type, value = VE,  metric = "VE")]
)
long_dt[, metric := factor(metric, levels = c("AVE", "VE"))]
long_dt[, model_type := factor(model_type, levels = models_to_run)]

# Reference lines (per-facet). Naive Wald applies to both facet-columns
# identically; posterior median/CI is per (metric, model_type).
ref_wald <- data.table(
  metric = factor(c("AVE", "AVE", "AVE", "VE", "VE", "VE"),
                  levels = c("AVE", "VE")),
  what   = c("center", "lo", "hi", "center", "lo", "hi"),
  value  = c(ave_naive, ave_naive_lo, ave_naive_hi,
             ve_naive,  ve_naive_lo,  ve_naive_hi))

ref_post <- rbindlist(lapply(names(diag_results), function(mt) {
  r <- diag_results[[mt]]
  rbind(
    data.table(model_type = mt, metric = "AVE", what = "center", value = r$ave_q[2]),
    data.table(model_type = mt, metric = "AVE", what = "lo",     value = r$ave_q[1]),
    data.table(model_type = mt, metric = "AVE", what = "hi",     value = r$ave_q[3]),
    data.table(model_type = mt, metric = "VE",  what = "center", value = r$ve_q[2]),
    data.table(model_type = mt, metric = "VE",  what = "lo",     value = r$ve_q[1]),
    data.table(model_type = mt, metric = "VE",  what = "hi",     value = r$ve_q[3]))
}))
ref_post[, metric := factor(metric, levels = c("AVE", "VE"))]
ref_post[, model_type := factor(model_type, levels = models_to_run)]

p <- ggplot(long_dt, aes(x = value)) +
  geom_histogram(bins = 40, fill = pal_post, alpha = 0.7) +
  # Wald reference (same in every facet within a row)
  geom_vline(data = ref_wald[what == "center"],
             aes(xintercept = value), colour = pal_naive, size = 0.9) +
  geom_vline(data = ref_wald[what != "center"],
             aes(xintercept = value), colour = pal_naive, linetype = "dashed") +
  # Posterior CI (per facet)
  geom_vline(data = ref_post[what == "center"],
             aes(xintercept = value), colour = pal_post, size = 0.9) +
  geom_vline(data = ref_post[what != "center"],
             aes(xintercept = value), colour = pal_post, linetype = "dashed") +
  facet_grid(metric ~ model_type, scales = "free") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = NULL, y = "count",
       title = glue("Posterior propagation (blue) vs Wald (orange) ",
                    "— study N_cont = {N_cont}, N_vac = {N_vac}, ",
                    "C1 = {data_C1}, C2 = {data_C2}, t* = {t_star}"))
ggsave(file.path(out_dir, "posterior_vs_wald.png"), p,
       width = 12, height = 8, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
