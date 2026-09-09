# Spread of C1, C2 and their ratio ACROSS ALLOCATIONS, at fixed parameters.
#
# pl_alpha = 1.4, alpha = 0.3, beta calibrated per network seed so that
# mean C1 ~ 160 of 400 unvaccinated at t_star = 8. Then hold (beta, alpha)
# fixed and vary only WHICH 400 of the 800 nodes are vaccinated.
#
# Each allocation's C1/C2 is a mean over n_sim realisations, so the observed
# across-allocation sd contains Monte-Carlo noise; we subtract it to get the
# true between-allocation component.

suppressMessages(library(data.table))
source("utils.R"); source("net_sir_events.R")
net_sir_compile()

N <- 800; N_vac <- 400; N_unvac <- N - N_vac
mean_k <- 6; pa <- 1.4
t_star <- 8; init_I <- 2; gamma <- 1
alpha_fix <- 0.3
target_C1 <- 160
n_alloc <- 30            # allocations per network seed
n_sim   <- 400           # realisations per allocation
net_seeds <- 1:5
cores <- 8

sim_alloc <- function(csr, vac, b, nsim, seed) {
  unvac <- setdiff(seq_len(N), vac)
  sus <- rep(1, N); sus[vac] <- alpha_fix
  inf <- run_stoch_network_events(
    beta = b, N = N, susceptibility = sus, t = t_star, vac = vac, csr = csr,
    gamma = gamma, timepoints = t_star, I_ini = init_I, n_sim = nsim,
    seed = seed, k_mean = mean_k, cores = cores, return_times = TRUE)
  hit <- inf <= t_star
  list(C1 = rowSums(hit[, unvac, drop = FALSE]),
       C2 = rowSums(hit[, vac,   drop = FALSE]))
}

res <- rbindlist(lapply(net_seeds, function(ns) {
  adj <- sample_pareto_adj(N, alpha = pa, mean_k = mean_k, seed = ns)
  csr <- adj_to_csr(adj = adj)
  iso <- mean(adj$degree == 0)

  # calibrate beta on the average over 5 reference allocations, so beta is
  # not tuned to any single one
  ref <- lapply(1:5, function(a) { set.seed(700 + a); sample.int(N, N_vac) })
  mC1 <- function(b) mean(vapply(seq_along(ref), function(a)
    mean(sim_alloc(csr, ref[[a]], b, 200, 3 + a)$C1), numeric(1)))
  lo <- log(0.2); hi <- log(60)
  for (i in 1:16) { m <- (lo + hi) / 2
    if (mC1(exp(m)) < target_C1) lo <- m else hi <- m }
  b <- exp((lo + hi) / 2)

  cat(sprintf("net seed %d: beta = %.3f  (isolated %.1f%%, mean deg %.2f)\n",
              ns, b, 100 * iso, mean(adj$degree)))
  flush.console()

  rbindlist(lapply(seq_len(n_alloc), function(a) {
    set.seed(2000 + 41 * ns + a); vac <- sample.int(N, N_vac)
    o <- sim_alloc(csr, vac, b, n_sim, 400 + a)
    tot <- o$C1 + o$C2
    data.table(net_seed = ns, beta = b, iso = iso, alloc = a,
               C1 = mean(o$C1), C2 = mean(o$C2),
               sd_C1 = sd(o$C1), sd_C2 = sd(o$C2),
               fizzle = mean(tot <= 0.05 * N),
               n_sim = n_sim)
  }))
}))

res[, `:=`(ratio_12 = C1 / C2, CIR = C2 / C1)]
fwrite(res, "output/alloc_spread.csv")

# between-allocation sd, with the Monte-Carlo component removed
decomp <- function(x, mc) {
  s <- sd(x); v <- s^2 - mean(mc^2)
  c(sd_obs = s, sd_mc = sqrt(mean(mc^2)), sd_true = sqrt(max(v, 0)))
}

cat(sprintf("\n=== pl_alpha=%.1f  alpha=%.2f  t*=%d  %d allocations x %d reps ===\n",
            pa, alpha_fix, t_star, n_alloc, n_sim))
for (ns in net_seeds) {
  d <- res[net_seed == ns]
  d1 <- decomp(d$C1, d$sd_C1 / sqrt(d$n_sim))
  d2 <- decomp(d$C2, d$sd_C2 / sqrt(d$n_sim))
  cat(sprintf("\n-- net seed %d (beta=%.2f, isolated %.1f%%, P(fizzle)=%.3f) --\n",
              ns, d$beta[1], 100 * d$iso[1], mean(d$fizzle)))
  cat(sprintf("  C1     mean %6.1f  sd %5.2f (MC %4.2f -> alloc %5.2f)  range %.0f-%.0f  cv %.3f\n",
              mean(d$C1), d1[["sd_obs"]], d1[["sd_mc"]], d1[["sd_true"]],
              min(d$C1), max(d$C1), d1[["sd_true"]] / mean(d$C1)))
  cat(sprintf("  C2     mean %6.1f  sd %5.2f (MC %4.2f -> alloc %5.2f)  range %.0f-%.0f  cv %.3f\n",
              mean(d$C2), d2[["sd_obs"]], d2[["sd_mc"]], d2[["sd_true"]],
              min(d$C2), max(d$C2), d2[["sd_true"]] / mean(d$C2)))
  cat(sprintf("  C1/C2  mean %6.3f  sd %.4f  range %.3f-%.3f  cv %.3f\n",
              mean(d$ratio_12), sd(d$ratio_12), min(d$ratio_12), max(d$ratio_12),
              sd(d$ratio_12) / mean(d$ratio_12)))
  cat(sprintf("  CIR    mean %6.3f  sd %.4f  range %.3f-%.3f  cv %.3f   (VE = %.3f)\n",
              mean(d$CIR), sd(d$CIR), min(d$CIR), max(d$CIR),
              sd(d$CIR) / mean(d$CIR), 1 - mean(d$CIR)))
}

cat("\n=== summary across network seeds ===\n")
print(res[, .(beta = round(beta[1], 2),
              C1 = round(mean(C1), 1), sd_C1 = round(sd(C1), 2),
              C2 = round(mean(C2), 1), sd_C2 = round(sd(C2), 2),
              CIR = round(mean(CIR), 4), sd_CIR = round(sd(CIR), 4),
              cv_CIR = round(sd(CIR) / mean(CIR), 4),
              fizzle = round(mean(fizzle), 3)), by = net_seed])

cat("\nNote: cv_CIR is the coefficient of variation of the trial contrast\n",
    "across allocations at FIXED parameters -- the raw allocation signal,\n",
    "before any compensating refit and before averaging over allocations.\n")
cat("\nWrote output/alloc_spread.csv\n")
