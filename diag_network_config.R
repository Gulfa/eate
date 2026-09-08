# Diagnose one network config: expA_N10000I20__network_pa1.2_n6_a2
# ---------------------------------------------------------------------------
# The fit itself looked fine (loss_chisq = 0.568), but the POSTERIOR did not:
#
#   grid MLE: beta = 2.921 alpha = 0.2366 | 95% LR CI beta [1.47, 5] alpha [0.18, 0.35]
#   grid diag: span=8.8 interior=FALSE bnd_mass=0.0115 ESS=78/625 sd/box=(0.59, 0.26)  [BOX-LIMITED]
#   sd_beta = 0.844
#
# The beta interval runs to 5, which is exp(log_beta_hi), i.e. the fitting
# bound -- so the likelihood is not closing in beta and the reported sd is
# partly the width of the box rather than the data. This script rebuilds that
# exact network and allocation and asks why.
#
#   Rscript diag_network_config.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table)})
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()

# ---- the config, from experiments[["expA_N10000I20"]] ---------------------
N_cont   <- 400; N_vac <- 400
N        <- N_cont + N_vac
data_C1  <- 160; data_C2 <- 80
t_star   <- 8
init_I   <- 2                 # init_I_nw
pl_alpha <- 1.2
mean_k   <- 6
gamma    <- 1
net_seed <- 6
alloc_seed <- 2
beta_fit <- 2.7162
alpha_fit <- 0.2389
n_rep    <- 1000
cores    <- 8

adj <- sample_pareto_adj(N, alpha = pl_alpha, mean_k = mean_k, seed = net_seed)
csr <- adj_to_csr(adj = adj)
set.seed(alloc_seed); vac <- sample.int(N, N_vac)
unvac <- setdiff(seq_len(N), vac)
tp <- seq(1, t_star, 1)

cat(sprintf("network: N=%d  pl_alpha=%.1f  seed=%d\n", N, pl_alpha, net_seed))
d <- adj$degree
cat(sprintf("  degree: mean=%.2f median=%d max=%d  q99=%.0f  isolated=%d (%.1f%%)\n",
            mean(d), median(d), max(d), quantile(d, .99), sum(d == 0),
            100 * mean(d == 0)))
cat(sprintf("  share of all edges held by the top 1%% of nodes: %.1f%%\n",
            100 * sum(sort(d, decreasing = TRUE)[1:round(N/100)]) / sum(d)))
cat(sprintf("  vaccinated: %d of %d; mean degree vac=%.2f unvac=%.2f\n\n",
            length(vac), N, mean(d[vac]), mean(d[unvac])))

sim <- function(b, a, nsim = n_rep, seed = 11) {
  sus <- rep(1, N); sus[vac] <- a
  inf <- run_stoch_network_events(
    beta = b, N = N, susceptibility = sus, t = t_star, vac = vac, csr = csr,
    gamma = gamma, timepoints = tp, I_ini = init_I, n_sim = nsim, seed = seed,
    k_mean = mean_k, cores = cores, return_times = TRUE)
  hit <- inf <= t_star
  data.table(C1 = rowSums(hit[, unvac, drop = FALSE]),
             C2 = rowSums(hit[, vac,   drop = FALSE]))
}

# ---- 1. outcome distribution at the fit -----------------------------------
o <- sim(beta_fit, alpha_fit)
cat(sprintf("=== %d reps at the fitted beta=%.4f alpha=%.4f (target %d / %d) ===\n",
            n_rep, beta_fit, alpha_fit, data_C1, data_C2))
cat(sprintf("  C1: mean=%6.1f sd=%6.1f  median=%6.1f  q10=%.0f q90=%.0f\n",
            mean(o$C1), sd(o$C1), median(o$C1), quantile(o$C1,.1), quantile(o$C1,.9)))
cat(sprintf("  C2: mean=%6.1f sd=%6.1f  median=%6.1f  q10=%.0f q90=%.0f\n",
            mean(o$C2), sd(o$C2), median(o$C2), quantile(o$C2,.1), quantile(o$C2,.9)))
fizz <- (o$C1 + o$C2) <= 0.05 * N
cat(sprintf("  P(fizzle, total <= 5%% of N) = %.3f   -> bimodal: %s\n",
            mean(fizz), if (mean(fizz) > 0.02 && mean(fizz) < 0.98) "YES" else "no"))
if (any(fizz) && any(!fizz))
  cat(sprintf("  takeoffs only: C1 mean=%.1f  C2 mean=%.1f  (n=%d)\n",
              mean(o$C1[!fizz]), mean(o$C2[!fizz]), sum(!fizz)))
cat("\n")

# ---- 2. how do the outcomes respond to beta? ------------------------------
# If mean C1 saturates, the likelihood cannot close in beta and the LR region
# runs to the fitting bound -- which is what BOX-LIMITED is reporting.
cat("=== response to beta (alpha fixed at the fit) ===\n")
cat(sprintf("%-8s %-9s %-9s %-9s %s\n", "beta", "mean C1", "mean C2", "P(fizzle)", "C2/C1"))
for (b in c(1.0, 1.5, 2.0, 2.7162, 4, 6, 10, 20)) {
  x <- sim(b, alpha_fit, nsim = 600, seed = 21)
  cat(sprintf("%-8.2f %-9.1f %-9.1f %-9.3f %.3f\n", b, mean(x$C1), mean(x$C2),
              mean((x$C1 + x$C2) <= 0.05 * N), mean(x$C2) / mean(x$C1)))
}
cat(sprintf("(fitting bound is beta <= 5; the LR CI ran to it)\n\n"))

# ---- 3. and to alpha? -----------------------------------------------------
cat("=== response to alpha (beta fixed at the fit) ===\n")
cat(sprintf("%-8s %-9s %-9s %s\n", "alpha", "mean C1", "mean C2", "C2/C1"))
for (a in c(0.05, 0.15, 0.2389, 0.35, 0.5, 0.8)) {
  x <- sim(beta_fit, a, nsim = 600, seed = 31)
  cat(sprintf("%-8.3f %-9.1f %-9.1f %.3f\n", a, mean(x$C1), mean(x$C2),
              mean(x$C2) / mean(x$C1)))
}

# ---- 4. is the beta bound the binding constraint? -------------------------
# Profile the SAME kernel likelihood the fit uses, out well past the current
# log_beta_hi (= log 5). If -log L is still falling at 5, the bound is
# truncating a region that would otherwise close, and BOX-LIMITED is simply
# saying "raise the bound". If it has already turned over well before 5, the
# region is wide because the data cannot pin beta, and a wider box will not
# help.
kbw <- c(max(0.5 * sqrt(data_C1 * (1 - data_C1 / N_cont)), 1),
         max(0.5 * sqrt(data_C2 * (1 - data_C2 / N_vac)),  1))
nll_at <- function(b, a, nsim = 4000, seed = 41) {
  x <- sim(b, a, nsim = nsim, seed = seed)
  k <- exp(-0.5 * (((x$C1 - data_C1) / kbw[1])^2 + ((x$C2 - data_C2) / kbw[2])^2))
  -log(mean(k) / (2 * pi * kbw[1] * kbw[2]) + 1e-300)
}
cat(sprintf("\n=== kernel -log L profile in beta (h = %.1f, %.1f), alpha at the fit ===\n",
            kbw[1], kbw[2]))
cat(sprintf("%-8s %-10s %-10s %s\n", "beta", "-log L", "P(fizzle)", ""))
bs <- c(1.5, 2, 2.7162, 3.5, 4, 5, 7, 10, 15, 25)
nl <- vapply(bs, function(b) nll_at(b, alpha_fit), numeric(1))
pf <- vapply(bs, function(b) { x <- sim(b, alpha_fit, 600, 51)
                               mean((x$C1 + x$C2) <= 0.05 * N) }, numeric(1))
for (i in seq_along(bs))
  cat(sprintf("%-8.2f %-10.2f %-10.3f %s\n", bs[i], nl[i], pf[i],
              if (abs(bs[i] - 5) < 1e-9) "<- current bound" else ""))
cat(sprintf("\nminimum of the profile at beta = %.2f\n", bs[which.min(nl)]))
cat("If -log L is still dropping at beta = 5, the BOUND is the problem.\n",
    "If it turned over earlier, the wide interval is the DATA, not the box.\n")
