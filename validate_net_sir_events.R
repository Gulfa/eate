# Validate the exact event-driven network SIR against the dust tau-leaping
# model, before letting it anywhere near the pipeline.
#
# The two are the same continuous-time process, so they must agree in
# DISTRIBUTION (not realisation-by-realisation -- the RNG streams differ).
# The dust version carries a dt discretisation bias, so agreement should
# IMPROVE as dt shrinks; that is the check that matters.
#
#   Rscript validate_net_sir_events.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(dplyr)})
source("utils.R")            # get_conact_matrix_pl
source("stoch_model.R")
source("net_sir_events.R")

N        <- 1000
k_mean   <- 6
pl_alpha <- 3
t_star   <- 8
gamma    <- 1
n_sim    <- 1500
dts      <- c(0.1, 0.02)                 # bias should shrink with dt
grid     <- list(c(beta = 1.5, alpha = 0.5),
                 c(beta = 2.5, alpha = 0.3),
                 c(beta = 4.0, alpha = 0.8))

set.seed(1); c_ij <- get_conact_matrix_pl(N, alpha = pl_alpha, mean_k = k_mean)
adj  <- contact_matrix_to_adj(c_ij)
csr  <- adj_to_csr(adj = adj)
set.seed(2); vac <- sample(seq_len(N), N / 2)
tp   <- seq(1, t_star, 1)

cat(sprintf("N=%d mean_deg=%.1f max_deg=%d edges=%d  n_sim=%d\n\n",
            N, mean(lengths(split(csr$nbr, rep(seq_len(N), diff(csr$ptr))))),
            adj$max_degree, length(csr$nbr), n_sim))

res <- rbindlist(lapply(grid, function(p) {
  # --- event-driven -------------------------------------------------------
  te <- system.time(ev <- run_stoch_network_events(
    beta = p[["beta"]], N = N, susceptibility = c(1, p[["alpha"]]),
    t = t_star, vac = vac, csr = csr, gamma = gamma, timepoints = tp,
    I_ini = 2, n_sim = n_sim, seed = 11, k_mean = k_mean))[["elapsed"]]
  e <- ev[time == t_star]

  rbindlist(lapply(dts, function(d) {
    td <- system.time(du <- run_stoch_network(
      beta = p[["beta"]], N = N, susceptibility = c(1, p[["alpha"]]),
      t = t_star, c_ij = c_ij, adj = adj, vac = vac, k_mean = k_mean,
      gamma = gamma, dt = d, timepoints = tp, I_ini = 2,
      n_sim = n_sim, cores = 4, seed = 7))[["elapsed"]]
    setDT(du); f <- du[time == t_star]
    ks1 <- suppressWarnings(ks.test(e$C1, f$C1)$p.value)
    ks2 <- suppressWarnings(ks.test(e$C2, f$C2)$p.value)
    data.table(beta = p[["beta"]], alpha = p[["alpha"]], dt = d,
               dust_C1 = mean(f$C1), ev_C1 = mean(e$C1),
               dust_C2 = mean(f$C2), ev_C2 = mean(e$C2),
               rel_C1 = mean(f$C1) / mean(e$C1) - 1,
               rel_C2 = mean(f$C2) / mean(e$C2) - 1,
               sd_C1_dust = sd(f$C1), sd_C1_ev = sd(e$C1),
               ks_C1 = ks1, ks_C2 = ks2,
               t_dust = td, t_event = te, speedup = td / te)
  }))
}))

cat("=== means (dust tau-leaping vs exact event-driven) ===\n")
print(res[, .(beta, alpha, dt,
              dust_C1 = round(dust_C1, 1), ev_C1 = round(ev_C1, 1),
              rel_C1  = sprintf("%+.1f%%", 100 * rel_C1),
              dust_C2 = round(dust_C2, 1), ev_C2 = round(ev_C2, 1),
              rel_C2  = sprintf("%+.1f%%", 100 * rel_C2))])

cat("\n=== distribution (KS p-value; small p = distributions differ) + spread ===\n")
print(res[, .(beta, alpha, dt,
              sd_C1_dust = round(sd_C1_dust, 1), sd_C1_ev = round(sd_C1_ev, 1),
              ks_C1 = signif(ks_C1, 3), ks_C2 = signif(ks_C2, 3))])

cat("\n=== timing ===\n")
print(res[, .(beta, alpha, dt, t_dust = round(t_dust, 2),
              t_event = round(t_event, 3), speedup = round(speedup, 1))])

cat("\nPASS if: rel_* shrinks towards 0 as dt shrinks (dust converging on the\n",
    "exact process), sd matches, and KS p-values are not tiny at the small dt.\n")
fwrite(res, "output/validate_net_sir_events.csv")
