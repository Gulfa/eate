# Does the per-allocation refit really remove the between-allocation spread
# in the ABSOLUTE effect (AVE)? Measured, not argued.
#
# For each allocation a_j of one network:
#   1. calibrate (beta_j, alpha_j) to hit C1=160, C2=80 CONDITIONAL on a_j
#      (this is what fit_one does for each job)
#   2. AVE_fresh / VE_fresh = evaluated on 5 FRESH allocations  <- pipeline
#   3. AVE_own   / VE_own   = evaluated on a_j itself
#   4. AVE_fixed / VE_fixed = evaluated on a_j at the MEAN parameters
#      (allocation effect with no compensating refit)
#
# If the refit is what kills the absolute spread, sd(AVE_fresh) << sd(AVE_fixed).
# If instead the refit re-expresses it, sd(AVE_fresh) should be comparable
# or larger.

suppressMessages(library(data.table))
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()

N <- 800; N_vac <- 400; mean_k <- 6; pa <- 1.4
t_star <- 8; init_I <- 2; gamma <- 1
data_C1 <- 160; data_C2 <- 80
J <- 12
net_seed <- 1

set.seed(net_seed)
cm  <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
adj <- contact_matrix_to_adj(cm)
csr <- adj_to_csr(adj = adj)
allocs <- lapply(seq_len(J), function(j) { set.seed(3000 + j); sample.int(N, N_vac) })

# --- fast arm outcomes via the event engine (same beta scale) --------------
arms <- function(b, a, vac, nsim, seed) {
  sus <- rep(1, N); sus[vac] <- a
  r <- run_stoch_network(beta = b, N = N, susceptibility = sus, t = t_star,
        vac = vac, gamma = gamma, timepoints = t_star, I_ini = init_I,
        n_sim = nsim, seed = seed, k_mean = mean_k, cores = 1,
        engine = "events", csr = csr)
  setDT(r); f <- r[time == t_star]
  c(C1 = mean(f$C1), C2 = mean(f$C2))
}

calibrate <- function(vac) {
  a <- 0.3; b <- 5
  for (round in 1:3) {
    lo <- log(0.2); hi <- log(60)
    for (i in 1:15) { m <- (lo + hi) / 2
      if (arms(exp(m), a, vac, 300, 11)[["C1"]] < data_C1) lo <- m else hi <- m }
    b <- exp((lo + hi) / 2)
    alo <- log(0.01); ahi <- log(2)
    for (i in 1:13) { m <- (alo + ahi) / 2
      if (arms(b, exp(m), vac, 300, 13)[["C2"]] < data_C2) alo <- m else ahi <- m }
    a <- exp((alo + ahi) / 2)
  }
  c(beta = b, alpha = a)
}

fits <- parallel::mclapply(allocs, calibrate, mc.cores = min(J, 12))
fit  <- as.data.table(do.call(rbind, fits))[, j := seq_len(J)]
cat("=== per-allocation fits (same network, same target data) ===\n")
print(fit[, .(j, beta = round(beta, 3), alpha = round(alpha, 4))])
cat(sprintf("\nbeta : mean %.3f sd %.3f cv %.3f   (range %.2f - %.2f)\n",
  mean(fit$beta), sd(fit$beta), sd(fit$beta)/mean(fit$beta), min(fit$beta), max(fit$beta)))
cat(sprintf("alpha: mean %.4f sd %.4f cv %.3f   (range %.3f - %.3f)\n\n",
  mean(fit$alpha), sd(fit$alpha), sd(fit$alpha)/mean(fit$alpha), min(fit$alpha), max(fit$alpha)))
flush.console()

b_bar <- mean(fit$beta); a_bar <- mean(fit$alpha)

# --- VE / AVE through the actual EATE function ----------------------------
eate_of <- function(b, a, vl, n_rep = 300) {
  r <- get_stoch_eate_network(beta = b, susceptibility = c(1, a), f = 0.5, N = N,
        t = t_star, c_ij = cm, adj = adj, k_mean = mean_k, gamma = gamma,
        n_rep = n_rep, timepoints = seq(1, t_star, 1), init_I = init_I,
        vac_list = vl, mc.cores = 1, inner_cores = 1)
  setDT(r)
  # full grid, then select t_star: .cum_trapz's first row is 0, so a scalar
  # timepoints would zero the counterfactual and collapse this to the arm ratio
  fs <- r[method == "full_stoch" & t == t_star]
  c(VE = mean(1 - fs$eate), AVE = mean(fs$ave))
}

out <- rbindlist(parallel::mclapply(seq_len(J), function(j) {
  vac <- allocs[[j]]
  fresh <- lapply(1:5, function(s) { set.seed(9000 + 31*j + s); sample.int(N, N_vac) })
  a_own   <- eate_of(fit$beta[j], fit$alpha[j], list(vac))
  a_fresh <- eate_of(fit$beta[j], fit$alpha[j], fresh)
  a_fixed <- eate_of(b_bar,       a_bar,        list(vac))
  chk <- arms(fit$beta[j], fit$alpha[j], vac, 600, 21)
  data.table(j = j, beta = fit$beta[j], alpha = fit$alpha[j],
    VE_own = a_own[["VE"]],   AVE_own = a_own[["AVE"]],
    VE_fresh = a_fresh[["VE"]], AVE_fresh = a_fresh[["AVE"]],
    VE_fixed = a_fixed[["VE"]], AVE_fixed = a_fixed[["AVE"]],
    C1_chk = chk[["C1"]], C2_chk = chk[["C2"]])
}, mc.cores = min(J, 12)))

cat("=== per-allocation VE / AVE ===\n")
print(out[, .(j, beta = round(beta,2), alpha = round(alpha,3),
              VE_own = round(VE_own,4), VE_fresh = round(VE_fresh,4),
              VE_fixed = round(VE_fixed,4),
              AVE_own = round(AVE_own,4), AVE_fresh = round(AVE_fresh,4),
              AVE_fixed = round(AVE_fixed,4),
              C1 = round(C1_chk), C2 = round(C2_chk))])

sh <- function(lab, x) cat(sprintf("  %-34s mean %7.4f  sd %7.4f  cv %6.3f  range %.3f - %.3f\n",
  lab, mean(x), sd(x), sd(x)/abs(mean(x)), min(x), max(x)))
cat("\n=== spread ACROSS ALLOCATIONS ===\n VE:\n")
sh("own allocation (fit+eval consistent)", out$VE_own)
sh("fresh allocations (pipeline)", out$VE_fresh)
sh("own allocation, params NOT refitted", out$VE_fixed)
cat(" AVE (absolute, per person):\n")
sh("own allocation (fit+eval consistent)", out$AVE_own)
sh("fresh allocations (pipeline)", out$AVE_fresh)
sh("own allocation, params NOT refitted", out$AVE_fixed)

cat(sprintf("\n  cor(beta_j, AVE_fresh) = %+.3f    cor(beta_j, AVE_own) = %+.3f\n",
            cor(out$beta, out$AVE_fresh), cor(out$beta, out$AVE_own)))
cat(sprintf("  sd(AVE_fresh)/sd(AVE_fixed) = %.2f   <- <1 means the refit REMOVES spread\n",
            sd(out$AVE_fresh)/sd(out$AVE_fixed)))
cat(sprintf("  sd(VE_fresh) /sd(VE_fixed)  = %.2f\n",
            sd(out$VE_fresh)/sd(out$VE_fixed)))
fwrite(out, "output/ave_refit_check.csv")
cat("\nWrote output/ave_refit_check.csv\n")
