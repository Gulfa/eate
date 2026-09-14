# How many index cases give a 50/50 split between fizzle and a 160/80 epidemic?
# ---------------------------------------------------------------------------
# For the bimodal regime the kernel fit was built for: half the realisations
# producing no epidemic, half producing roughly the target case counts.
#
# The two knobs turn out to be independent. Conditional on takeoff the final
# size is set by beta and alpha, not by the seed count -- across I_ini = 1..12
# the takeoff-conditional counts move only 157-163 / 79-82 -- so calibrate
# beta/alpha once on the takeoff-conditional means and then tune I_ini purely
# for the fizzle probability.
#
# CAVEAT: run_stoch_network_events seeds the FIRST I_ini nodes deterministically
# (net_sir_events.R, `seeds <- seq_len(min(I_ini, n)) - 1L`), so the fizzle
# probability depends on whether those particular nodes happen to be hubs, and
# varies a lot between network realisations -- 0.44 to 0.69 across four networks
# at I_ini = 3, and at I_ini = 1 two of the four never take off at all. Hence
# the scan over net_seeds. Tune per network seed if you need a specific one at
# 50%.
#
#   Rscript diag_seed_fizzle.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({library(data.table)}))
source("utils.R"); suppressWarnings(suppressMessages(source("stoch_model.R")))
source("net_sir_events.R"); net_sir_compile()

N <- 800; N_vac <- 400; mean_k <- 6; pa <- 1.5
t_star <- 8; gamma <- 1
data_C1 <- 160; data_C2 <- 80
fizz_cut <- 0.05 * N                 # total cases below this = no epidemic
nsim <- 2000; ref_I <- 3             # calibrate at this, then scan
inits <- c(1, 2, 3, 4, 6, 8, 12)
net_seeds <- 1:4

sim <- function(csr, vac, unvac, b, a, I0, nsim. = nsim, seed = 7) {
  sus <- rep(1, N); sus[vac] <- a
  inf <- run_stoch_network_events(beta = b, N = N, susceptibility = sus,
          t = t_star, vac = vac, csr = csr, gamma = gamma, timepoints = t_star,
          I_ini = I0, n_sim = nsim., seed = seed, k_mean = mean_k, cores = 8,
          return_times = TRUE)
  hit <- inf <= t_star
  list(C1 = rowSums(hit[, unvac, drop = FALSE]),
       C2 = rowSums(hit[, vac,   drop = FALSE]))
}

# calibrate on the TAKEOFF-conditional means: the target 160/80 describes an
# epidemic that happened, not an average over runs that mostly fizzled
cal <- function(csr, vac, unvac) {
  a <- 0.3
  for (rd in 1:3) {
    lo <- log(0.3); hi <- log(60)
    for (i in 1:14) { m <- (lo + hi) / 2
      o <- sim(csr, vac, unvac, exp(m), a, ref_I, 600, 3)
      k <- (o$C1 + o$C2) > fizz_cut
      v <- if (sum(k) > 5) mean(o$C1[k]) else 0
      if (v < data_C1) lo <- m else hi <- m }
    b <- exp((lo + hi) / 2)
    alo <- log(0.01); ahi <- log(2)
    for (i in 1:12) { m <- (alo + ahi) / 2
      o <- sim(csr, vac, unvac, b, exp(m), ref_I, 600, 3)
      k <- (o$C1 + o$C2) > fizz_cut
      v <- if (sum(k) > 5) mean(o$C2[k]) else 0
      if (v < data_C2) alo <- m else ahi <- m }
    a <- exp((alo + ahi) / 2)
  }
  c(beta = b, alpha = a)
}

res <- rbindlist(lapply(net_seeds, function(ns) {
  set.seed(ns); cm <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
  adj <- contact_matrix_to_adj(cm); csr <- adj_to_csr(contact_matrix = cm, adj = adj)
  set.seed(100 + ns); vac <- sample.int(N, N_vac)
  unvac <- setdiff(seq_len(N), vac)
  p <- cal(csr, vac, unvac)
  rbindlist(lapply(inits, function(I0) {
    o <- sim(csr, vac, unvac, p[["beta"]], p[["alpha"]], I0, nsim, 11)
    tot <- o$C1 + o$C2; k <- tot > fizz_cut
    data.table(net_seed = ns, beta = p[["beta"]], alpha = p[["alpha"]], I_ini = I0,
               p_fizzle = mean(!k),
               C1_takeoff = if (any(k)) mean(o$C1[k]) else NA_real_,
               C2_takeoff = if (any(k)) mean(o$C2[k]) else NA_real_,
               deg_seeds = sum(adj$degree[seq_len(I0)]))
  }))
}))
fwrite(res, "output/seed_scan_pa15.csv")

cat(sprintf("pa = %.1f, N = %d, t* = %d, target %d/%d conditional on takeoff\n",
            pa, N, t_star, data_C1, data_C2))
cat(sprintf("fizzle = total cases <= %d (5%% of N); %d sims per cell\n\n", fizz_cut, nsim))
cat("fitted per network:\n")
print(unique(res[, .(net_seed, beta = round(beta, 2), alpha = round(alpha, 3))]))
cat("\nP(no epidemic) by number of index cases:\n")
print(dcast(res, I_ini ~ net_seed, value.var = "p_fizzle")[,
  lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("\nmean over networks, with the takeoff-conditional case counts:\n")
print(res[, .(p_fizzle = round(mean(p_fizzle), 3),
              C1_takeoff = round(mean(C1_takeoff, na.rm = TRUE)),
              C2_takeoff = round(mean(C2_takeoff, na.rm = TRUE))), by = I_ini])
cat("\nWrote output/seed_scan_pa15.csv\n")
