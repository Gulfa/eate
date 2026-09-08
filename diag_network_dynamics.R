# Do the DYNAMICS depend on mean_k, or only the feasibility?
# ---------------------------------------------------------------------------
# diag_network_isolation.R showed that raising mean_k from 6 to 12 removes the
# isolated-node problem that made some heavy-tailed configs unfittable. But
# mean_k is not a free knob: it changes the contact structure, so it can change
# the epidemiology and hence the reported effects.
#
# Here each (pl_alpha, mean_k) cell is CALIBRATED to the same observed data --
# beta and alpha chosen so the mean (C1, C2) hits the target -- so every cell
# fits identically well and any remaining difference is the network structure
# talking. Reported per cell:
#
#   beta, alpha  the parameters needed to reproduce the same data
#   P(fizzle)    how often the epidemic dies out
#   averted/1k   infections averted by +10% coverage (the policy quantity,
#                panel B of the summary figure)
#
#   Rscript diag_network_dynamics.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()

N_cont <- 400; N_vac <- 400; N <- N_cont + N_vac
data_C1 <- 160; data_C2 <- 80
t_star <- 8; init_I <- 2; gamma <- 1
n_sim  <- 600
cores  <- 8
seeds  <- 1:4
pl_alphas <- c(1.2, 1.5, 2, 3)
mean_ks   <- c(6, 12)

sim_cov <- function(csr, vac, unvac, b, a, mk, nsim = n_sim, seed = 5) {
  sus <- rep(1, N); sus[vac] <- a
  inf <- run_stoch_network_events(
    beta = b, N = N, susceptibility = sus, t = t_star, vac = vac, csr = csr,
    gamma = gamma, timepoints = t_star, I_ini = init_I, n_sim = nsim,
    seed = seed, k_mean = mk, cores = cores, return_times = TRUE)
  hit <- inf <= t_star
  list(C1 = rowSums(hit[, unvac, drop = FALSE]),
       C2 = rowSums(hit[, vac,   drop = FALSE]))
}

# Alternate bisections: beta to match C1, alpha to match the C2/C1 ratio.
calibrate <- function(csr, vac, unvac, mk) {
  a <- 0.3
  for (round in 1:3) {
    lo <- log(0.05); hi <- log(200)
    for (i in 1:16) { m <- (lo + hi) / 2
      if (mean(sim_cov(csr, vac, unvac, exp(m), a, mk, 300, 3)$C1) < data_C1)
        lo <- m else hi <- m }
    b <- exp((lo + hi) / 2)
    alo <- log(0.005); ahi <- log(3)
    for (i in 1:14) { m <- (alo + ahi) / 2
      if (mean(sim_cov(csr, vac, unvac, b, exp(m), mk, 300, 3)$C2) < data_C2)
        alo <- m else ahi <- m }
    a <- exp((alo + ahi) / 2)
  }
  c(beta = b, alpha = a)
}

res <- rbindlist(lapply(mean_ks, function(mk)
  rbindlist(lapply(pl_alphas, function(pa)
    rbindlist(lapply(seeds, function(sd_) {
      adj <- sample_pareto_adj(N, alpha = pa, mean_k = mk, seed = sd_)
      csr <- adj_to_csr(adj = adj)
      set.seed(100 + sd_); vac <- sample.int(N, N_vac)
      unvac <- setdiff(seq_len(N), vac)
      iso <- mean(adj$degree == 0)

      p <- calibrate(csr, vac, unvac, mk)
      o <- sim_cov(csr, vac, unvac, p[["beta"]], p[["alpha"]], mk, n_sim, 7)
      fz <- mean((o$C1 + o$C2) <= 0.05 * N)

      # +10% coverage, same parameters, fresh allocations at each level
      inf_at <- function(cv) {
        mean(vapply(1:3, function(a2) {
          set.seed(500 + 17 * a2); v <- sample.int(N, round(cv * N))
          x <- sim_cov(csr, v, setdiff(seq_len(N), v), p[["beta"]], p[["alpha"]],
                       mk, 400, 9 + a2)
          mean(x$C1 + x$C2) }, numeric(1)))
      }
      av <- (inf_at(0.5) - inf_at(0.6)) / N * 1000

      data.table(mean_k = mk, pl_alpha = pa, seed = sd_, iso = iso,
                 beta = p[["beta"]], alpha = p[["alpha"]],
                 C1 = mean(o$C1), C2 = mean(o$C2), fizzle = fz,
                 averted_per1k = av)
    }))))))

fwrite(res, "output/network_dynamics.csv")

summ <- res[, .(iso = round(median(iso), 3),
                beta = round(median(beta), 2), alpha = round(median(alpha), 3),
                C1 = round(median(C1)), C2 = round(median(C2)),
                fizzle = round(median(fizzle), 3),
                averted = round(median(averted_per1k), 1),
                averted_min = round(min(averted_per1k), 1),
                averted_max = round(max(averted_per1k), 1)),
            by = .(mean_k, pl_alpha)][order(pl_alpha, mean_k)]
cat(sprintf("=== calibrated to C1=%d C2=%d, %d seeds per cell ===\n",
            data_C1, data_C2, length(seeds)))
print(summ)

res[, mk := factor(mean_k, labels = paste0("mean_k = ", mean_ks))]
long <- melt(res, id.vars = c("mk", "pl_alpha", "seed"),
             measure.vars = c("beta", "alpha", "fizzle", "averted_per1k"),
             variable.name = "quantity")
long[, quantity := factor(quantity,
        levels = c("beta", "alpha", "fizzle", "averted_per1k"),
        labels = c("fitted beta", "fitted alpha", "P(fizzle)",
                   "infections averted / 1000 (+10% coverage)"))]
p <- ggplot(long, aes(factor(pl_alpha), value, colour = mk)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.7),
             alpha = 0.7, size = 1.6) +
  facet_wrap(~ quantity, scales = "free_y") +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = "Pareto exponent (pl_alpha)", y = NULL,
       title = "Network structure changes the dynamics, not just the feasibility",
       subtitle = sprintf("every cell calibrated to the same data (C1 = %d, C2 = %d), so differences are structural",
                          data_C1, data_C2))
ggsave("output/network_dynamics.png", p, width = 10, height = 7, dpi = 140)
cat("\nWrote output/network_dynamics.{csv,png}\n")
