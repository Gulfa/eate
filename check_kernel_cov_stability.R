# check_kernel_cov_stability.R
# ---------------------------------------------------------------------------
# Is the kernel-Hessian posterior stable, or is it estimator noise?
#
# For each SIR I_ini setting: fit ONCE (kernel/mode), then recompute
# kernel_posterior_cov at that SAME (beta_hat, alpha_hat) under several
# post_cov_seed values. The fit is held fixed, so any spread across seeds is
# pure Monte-Carlo noise in the Hessian estimator.
#
#   stable across seeds  -> Sigma is real; non-monotonic sd_VE is a genuine
#                           quadratic-form effect (gradient + cross term).
#   jumps across seeds   -> estimator noise; raise kernel_hess_nsim (only the
#                           ~p_takeoff fraction of particles carries signal).
#
# Also reports the delta-method decomposition of Var(VE):
#   Var(VE) ~ g_b^2 s_b^2 + 2 g_b g_a s_ba + g_a^2 s_a^2,  g = dVE/d(beta,alpha)
# so you can see WHICH term drives the non-monotonicity.
#
# Run (fhi container, from /fhi/eate):  Rscript check_kernel_cov_stability.R
# ---------------------------------------------------------------------------

Sys.setenv(EATE_SOURCE_ONLY = "1")
suppressWarnings(suppressMessages(source("run_fit_array.R")))
suppressMessages(library(data.table))

# ---- config ---------------------------------------------------------------
I_ini_list <- list(c(1, 1), c(5, 5), c(10, 10))   # sweep to check
cov_seeds  <- c(1234L, 5678L, 9012L)              # seeds for the stability check
cores      <- 8
exp_idx    <- 1                                   # which experiment's sizes/data

exp0 <- experiments[[exp_idx]]
base <- modifyList(base_common, list(
  experiment_id = exp0$id,
  N_cont = exp0$N_cont, N_vac = exp0$N_vac,
  data_C1 = exp0$data_C1, data_C2 = exp0$data_C2,
  t_star = exp0$t_star, I_ini_2g = exp0$I_ini_2g, init_I_nw = exp0$init_I_nw,
  model_type = "sir", ve_n_vac = 1, inner_cores = cores,
  fit_method = "kernel"))

# VE(beta, alpha) at t*, averaged over inner sims (for the delta-method gradient).
ve_at <- function(cfg, b, a) {
  v <- compute_ve(cfg, b, a); setDT(v)
  mean(1 - v[method == "full_stoch" & t == cfg$t_star, eate], na.rm = TRUE)
}

res <- rbindlist(lapply(I_ini_list, function(iv) {
  cfg <- modifyList(base, list(I_ini_2g = iv))
  sim <- build_simulator(cfg)
  fit <- fit_one(sim, cfg)

  # takeoff fraction at the fit (how many particles actually inform the kernel)
  o  <- sim(fit$beta, fit$alpha, 2000, seed = 7)
  pt <- mean((o$C1 + o$C2) > 0.05 * (cfg$N_cont + cfg$N_vac))

  # same (beta_hat, alpha_hat), several cov seeds -> pure estimator noise
  rows <- rbindlist(lapply(cov_seeds, function(s) {
    cfg_s <- modifyList(cfg, list(post_cov_seed = s))
    kc <- kernel_posterior_cov(sim, cfg_s, fit$beta, fit$alpha)
    cv <- kc$cov
    data.table(seed = s, sd_b = kc$sd[["beta"]], sd_a = kc$sd[["alpha"]],
               rho = cv[1, 2] / sqrt(cv[1, 1] * cv[2, 2]))
  }))

  # delta-method gradient of VE at the fit (central differences)
  db <- 0.05 * fit$beta; da <- 0.05 * fit$alpha
  g_b <- (ve_at(cfg, fit$beta + db, fit$alpha) - ve_at(cfg, fit$beta - db, fit$alpha)) / (2 * db)
  g_a <- (ve_at(cfg, fit$beta, fit$alpha + da) - ve_at(cfg, fit$beta, fit$alpha - da)) / (2 * da)

  # variance decomposition using the MEDIAN cov across seeds
  sb <- median(rows$sd_b); sa <- median(rows$sd_a); rh <- median(rows$rho)
  t_bb <- g_b^2 * sb^2
  t_ba <- 2 * g_b * g_a * rh * sb * sa
  t_aa <- g_a^2 * sa^2
  data.table(
    I_ini    = sum(iv),
    beta     = round(fit$beta, 3),
    alpha    = round(fit$alpha, 3),
    p_takeoff= round(pt, 3),
    sd_b_med = round(sb, 4),
    sd_b_cv  = round(sd(rows$sd_b) / mean(rows$sd_b), 3),   # across-seed noise
    sd_a_med = round(sa, 4),
    sd_a_cv  = round(sd(rows$sd_a) / mean(rows$sd_a), 3),
    rho_med  = round(rh, 3),
    rho_rng  = round(diff(range(rows$rho)), 3),
    g_b = round(g_b, 4), g_a = round(g_a, 4),
    var_bb = signif(t_bb, 3), var_cross = signif(t_ba, 3), var_aa = signif(t_aa, 3),
    sd_VE_delta = round(sqrt(pmax(t_bb + t_ba + t_aa, 0)), 4))
}))

cat("\n=== Kernel posterior stability across cov seeds, per I_ini ===\n")
print(res[, .(I_ini, beta, alpha, p_takeoff,
              sd_b_med, sd_b_cv, sd_a_med, sd_a_cv, rho_med, rho_rng)])
cat("  *_cv = across-seed coefficient of variation at a FIXED fit.\n")
cat("  cv << 0.1 and small rho_rng => estimator stable (Sigma is real).\n")
cat("  cv large / rho_rng large     => Hessian noise; raise kernel_hess_nsim.\n")

cat("\n=== Delta-method decomposition of Var(VE) ===\n")
print(res[, .(I_ini, g_b, g_a, var_bb, var_cross, var_aa, sd_VE_delta)])
cat("  Var(VE) ~ var_bb + var_cross + var_aa; var_cross is usually NEGATIVE\n")
cat("  (beta-alpha ridge) -- that cancellation is why sd_VE need not track sd_beta.\n")

dir.create("output", showWarnings = FALSE)
fwrite(res, "output/kernel_cov_stability.csv")
cat("\nWrote output/kernel_cov_stability.csv\n")
