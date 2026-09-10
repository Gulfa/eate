# Does the EATE VE actually move with coverage, as saturation says it must?
# ---------------------------------------------------------------------------
# With a frozen cumulative force of infection Lambda,
#     P_unvac = 1 - exp(-Lambda)      P_vac = 1 - exp(-alpha*Lambda)
#     VE      = 1 - P_vac / P_unvac
# so VE -> 1 - alpha as Lambda -> 0 and VE -> 0 as Lambda -> infinity. Coverage
# moves Lambda (more vaccinated -> smaller epidemic -> smaller Lambda), so VE
# must rise with coverage, steeply wherever the attack rate is high.
#
# This calibrates the SIR to the usual target (AR_unvac = 0.40, AR_vac = 0.20 at
# 50% coverage, i.e. C1 = 160 / C2 = 80 of 400) and then sweeps coverage,
# reporting the realised attack rate, the EATE VE, and the VE the closed form
# above predicts from the realised Lambda. If the code is right the two agree.
#
#   Rscript diag_ve_vs_coverage.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R")

N <- 800; t_star <- 8; gamma <- 1; dt <- 0.05
I_ini_2g <- c(3, 3); target_AR_unvac <- 0.40; target_AR_vac <- 0.20
n_rep <- 400; n_vac_alloc <- 3; cores <- 8
covs <- c(0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95)

# arm attack rates at a given coverage, with I_ini split PROPORTIONALLY to the
# group sizes -- a fixed c(3,3) would mean something different at every coverage
# (3 seeds in a group of 40 at 95% coverage), which would confound the sweep.
arms <- function(b, a, f, nsim = n_rep, seed = 3, prop_seed = TRUE) {
  n_u <- round(N * (1 - f)); n_v <- N - n_u
  ii  <- if (prop_seed) {
    tot <- sum(I_ini_2g); c(round(tot * (1 - f)), tot - round(tot * (1 - f)))
  } else I_ini_2g
  ii <- pmin(ii, c(n_u, n_v))
  r <- run_stoch_cd_dust(matrix(1, 2, 2), beta = b, N = c(n_u, n_v), t = t_star,
        I_ini = ii, susceptibility = c(1, a), gamma = gamma, dt = dt,
        timepoints = t_star, n_sim = nsim, cores = cores, seed = seed)
  setDT(r); fin <- r[time == t_star]
  c(AR_unvac = mean(fin$C1) / n_u, AR_vac = mean(fin$C2) / n_v)
}

ve_at <- function(b, a, f) {
  n_u <- round(N * (1 - f)); tot <- sum(I_ini_2g)
  ii  <- c(round(tot * (1 - f)), tot - round(tot * (1 - f)))
  ii  <- pmin(ii, c(n_u, N - n_u))
  # timepoints MUST be the full grid, not the scalar t_star: the counterfactual
  # comes from .cum_trapz(), whose first row is 0 by construction, so a
  # single-element grid makes the cumulative FOI identically zero and the EATE
  # silently collapses to the factual arm ratio.
  r <- get_stoch_eate_sir(beta = b, susceptibility = c(1, a), f = f, N = N,
        t = t_star, gamma = gamma, I_ini = ii, n_vac = n_vac_alloc,
        n_rep = n_rep, dt = dt, timepoints = seq(1, t_star, 1), mc.cores = cores)
  setDT(r); fs <- r[method == "full_stoch" & t == t_star]
  c(VE = mean(1 - fs$eate), AVE = mean(fs$ave))
}

# ---- calibrate at 50% coverage -------------------------------------------
a <- 0.4; b <- 2
for (round_ in 1:3) {
  lo <- log(0.2); hi <- log(20)
  for (i in 1:16) { m <- (lo + hi) / 2
    if (arms(exp(m), a, 0.5)[["AR_unvac"]] < target_AR_unvac) lo <- m else hi <- m }
  b <- exp((lo + hi) / 2)
  alo <- log(0.01); ahi <- log(1.5)
  for (i in 1:14) { m <- (alo + ahi) / 2
    if (arms(b, exp(m), 0.5)[["AR_vac"]] < target_AR_vac) alo <- m else ahi <- m }
  a <- exp((alo + ahi) / 2)
}
cat(sprintf("calibrated at 50%% coverage: beta = %.4f  alpha = %.4f\n", b, a))
chk <- arms(b, a, 0.5, 800, 5)
cat(sprintf("  check: AR_unvac = %.3f (target %.2f)  AR_vac = %.3f (target %.2f)\n\n",
            chk[["AR_unvac"]], target_AR_unvac, chk[["AR_vac"]], target_AR_vac))

res <- rbindlist(lapply(covs, function(f) {
  ar <- arms(b, a, f, 800, 7)
  ve <- ve_at(b, a, f)
  # Lambda implied by the realised unvaccinated attack rate, and the closed-form
  # VE it predicts -- the independent check on the simulated EATE.
  lam <- -log(1 - ar[["AR_unvac"]])
  ve_theory <- 1 - (1 - exp(-a * lam)) / (1 - exp(-lam))
  data.table(coverage = f, AR_unvac = ar[["AR_unvac"]], AR_vac = ar[["AR_vac"]],
             Lambda = lam, VE = ve[["VE"]], VE_closed_form = ve_theory,
             AVE = ve[["AVE"]])
}))

cat("=== VE vs coverage (beta, alpha fixed at the 50% calibration) ===\n")
print(res[, .(coverage, AR_unvac = round(AR_unvac, 3), AR_vac = round(AR_vac, 3),
              Lambda = round(Lambda, 3), VE = round(VE, 4),
              VE_closed = round(VE_closed_form, 4), AVE = round(AVE, 4))])
cat(sprintf("\nVE range across coverage: %.4f -> %.4f  (span %.4f)\n",
            min(res$VE), max(res$VE), max(res$VE) - min(res$VE)))
cat(sprintf("1 - alpha (the Lambda -> 0 limit) = %.4f\n", 1 - a))
cat(sprintf("max |VE - closed form| = %.4f\n",
            max(abs(res$VE - res$VE_closed_form))))

fwrite(res, "output/ve_vs_coverage.csv")
long <- melt(res, id.vars = "coverage",
             measure.vars = c("VE", "VE_closed_form", "AR_unvac", "AVE"))
p <- ggplot(res, aes(coverage)) +
  geom_line(aes(y = VE, colour = "EATE VE (simulated)"), linewidth = 1) +
  geom_point(aes(y = VE, colour = "EATE VE (simulated)"), size = 2) +
  geom_line(aes(y = VE_closed_form, colour = "closed form from Lambda"),
            linewidth = 0.8, linetype = "dashed") +
  geom_line(aes(y = AR_unvac, colour = "attack rate, unvaccinated"),
            linewidth = 0.8) +
  geom_hline(yintercept = 1 - a, linetype = "dotted", colour = "grey50") +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  scale_x_continuous(labels = scales::percent) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(x = "vaccine coverage", y = NULL,
       title = "VE rises with coverage because coverage lowers the force of infection",
       subtitle = sprintf("dotted = 1 - alpha = %.2f, the Lambda -> 0 ceiling", 1 - a))
ggsave("output/ve_vs_coverage.png", p, width = 9, height = 6, dpi = 140)
cat("\nWrote output/ve_vs_coverage.{csv,png}\n")
