# Simple RCT, fixed vaccine effect: does CIR track the causal VE across coverage?
# ---------------------------------------------------------------------------
# N people on a contact network, a fraction `coverage` vaccinated at random
# with a FIXED alpha (no contact-dependence). At each coverage:
#
#   VE_trial = 1 - CIR = 1 - AR_vac / AR_unvac      what the trial observes
#   VE_flip  = 1 - EATE                             causal, by re-simulating
#                                                   each individual's flip
#
# Coverage is swept 5% -> 95%. beta is CALIBRATED PER COVERAGE so the
# unvaccinated attack rate is ~target_ar at every level: otherwise raising
# coverage shrinks the epidemic (AR fell 0.53 -> 0.05 in an earlier run), and
# VE moves with saturation rather than with the mechanism -- and at high
# coverage there is barely an epidemic left to measure.
#
# The vaccine effect is CONTACT-DEPENDENT. Two shapes, set by `thresh`:
#   thresh = 0  smooth ramp:   alpha_eff = 1 - f * (1 - alpha)
#   thresh > 0  hard cutoff:   alpha_eff = alpha if f >= thresh, else 1
# with f[i] the fraction vaccinated in i's local neighbourhood (i included).
#
# The cutoff is the interesting case: a step means flipping ONE person can push
# a neighbour ACROSS the threshold, discontinuously switching that neighbour's
# protection on. The trial contrast compares vaccinated with unvaccinated and
# cannot see that spillover, so CIR and the causal flip-VE should come apart --
# unlike the smooth ramp, where they agree at every coverage.
#
#   Rscript rct_cir_vs_ve.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()

# ---- config ---------------------------------------------------------------
N         <- 1000
pl_alpha  <- 3
mean_k    <- 6
alpha_vac <- 0.1          # full effect VE = 0.9 at f = 1
beta      <- 2.0          # starting point for the per-coverage calibration
target_ar <- 0.5          # unvaccinated attack rate to hold fixed
gamma     <- 1
t_star    <- 15
init_I    <- 5
n_rep     <- 2000         # replicates for the factual arm
n_rep_ve  <- 400          # replicates per re-simulated counterfactual
n_flip    <- 150          # individuals flipped per coverage
cores     <- 8
coverages <- seq(0.05, 0.95, by = 0.10)
n_reps    <- 5            # independent trials (fresh allocation + RNG) per coverage
thresh    <- 0.5          # 0 = smooth ramp; >0 = hard cutoff at this local coverage

set.seed(1)
adj <- sample_pareto_adj(N, alpha = pl_alpha, mean_k = mean_k, seed = 1)
csr <- adj_to_csr(adj = adj)
tp  <- seq(1, t_star, 1)
cat(sprintf("N=%d mean_deg=%.1f alpha=%.2f beta=%.1f t*=%d  effect=%s\n\n",
            N, mean(adj$degree), alpha_vac, beta, t_star,
            if (thresh > 0) sprintf("hard cutoff at f=%.2f", thresh) else "smooth ramp"))

run_cov <- function(cov, rep = 1L) {
  set.seed(2000 + round(1000 * cov) + 100000L * rep)
  n_v   <- max(1L, round(cov * N))
  vac   <- sample.int(N, n_v)
  unvac <- setdiff(seq_len(N), vac)

  # Contact-dependent effect: alpha_eff depends on local coverage.
  sus <- vacfrac_susceptibility(vac, alpha_vac, adj = adj,
                                vac_frac_ref = 1, vac_frac_power = 1,
                                vac_frac_thresh = thresh)

  ar_unvac_at <- function(b, nsim, sd) {
    inf <- run_stoch_network_events(
      beta = b, N = N, susceptibility = sus, t = t_star, vac = vac, csr = csr,
      gamma = gamma, timepoints = tp, I_ini = init_I, n_sim = nsim,
      seed = sd + 1000L * rep,
      k_mean = mean_k, cores = cores, return_times = TRUE)
    mean(rowMeans((inf <= t_star)[, unvac, drop = FALSE]))
  }
  # Calibrate beta so the CONTROL arm sees ~target_ar at every coverage.
  # Bisection on log(beta): AR is monotone in beta and cheap on this engine.
  lo <- log(0.2); hi <- log(50)
  for (it in 1:14) {
    mid <- (lo + hi) / 2
    if (ar_unvac_at(exp(mid), 400, 3) < target_ar) lo <- mid else hi <- mid
  }
  beta_c <- exp((lo + hi) / 2)

  inf <- run_stoch_network_events(
    beta = beta_c, N = N, susceptibility = sus, t = t_star, vac = vac, csr = csr,
    gamma = gamma, timepoints = tp, I_ini = init_I, n_sim = n_rep,
    seed = 7L + 1000L * rep,
    k_mean = mean_k, cores = cores, return_times = TRUE)
  hit  <- inf <= t_star
  ar_v <- mean(rowMeans(hit[, vac,   drop = FALSE]))
  ar_u <- mean(rowMeans(hit[, unvac, drop = FALSE]))

  # Causal flip-VE: re-simulate each sampled individual with status flipped.
  ve <- get_stoch_eate_network_vacfrac(
    beta = beta_c, alpha = alpha_vac, N = N, t = t_star, adj = adj, csr = csr,
    vac_frac_power = 1, vac_frac_ref = 1, vac_frac_thresh = thresh,
    k_mean = mean_k, gamma = gamma,
    n_rep = n_rep_ve, n_flip = n_flip, timepoints = tp, init_I = init_I,
    mc.cores = 1, inner_cores = cores, vac_list = list(vac),
    engine = "events", crn_seed = 13L + 100L * rep)
  setDT(ve)

  vi <- integer(N); vi[vac] <- 1L
  nv <- colSums(adj$mask * matrix(vi[adj$neighbors], nrow(adj$neighbors), N))
  f_bar <- mean(((vi + nv) / (colSums(adj$mask) + 1))[vac])

  data.table(coverage = cov, rep = rep, beta = beta_c, f_bar = f_bar,
             mean_alpha_eff = mean(sus[vac]),
             AR_vac = ar_v, AR_unvac = ar_u,
             CIR = ar_v / ar_u, VE_trial = 1 - ar_v / ar_u,
             VE_flip = 1 - ve[method == "full_stoch" & t == t_star, eate])
}

raw <- rbindlist(lapply(coverages, function(cv)
  rbindlist(lapply(seq_len(n_reps), function(rp) run_cov(cv, rp)))))
raw[, gap := VE_flip - VE_trial]
fwrite(raw, "output/rct_cir_vs_ve_reps.csv")

# Average over reps; se_gap is the standard error of the mean gap, i.e. the
# yardstick for whether a gap is real or just allocation/simulation noise.
res <- raw[, .(beta = mean(beta), f_bar = mean(f_bar),
               mean_alpha_eff = mean(mean_alpha_eff),
               AR_unvac = mean(AR_unvac), AR_vac = mean(AR_vac),
               CIR = mean(CIR),
               VE_trial = mean(VE_trial), VE_flip = mean(VE_flip),
               gap = mean(gap), se_gap = sd(gap) / sqrt(.N)), by = coverage]

cat("=== trial contrast vs causal flip-VE, by coverage ===\n")
print(res[, .(coverage, beta = round(beta, 2), f_bar = round(f_bar, 3),
              a_eff    = round(mean_alpha_eff, 3),
              AR_unvac = round(AR_unvac, 3), AR_vac = round(AR_vac, 3),
              CIR      = round(CIR, 3),
              VE_trial = round(VE_trial, 3),
              VE_flip  = round(VE_flip, 3),
              gap      = round(gap, 3),
              se_gap   = round(se_gap, 3),
              t        = round(gap / se_gap, 1))])
cat(sprintf("\nmean gap = %+.4f  (pooled se %.4f, %d reps x %d coverages)\n",
            mean(res$gap), sd(raw$gap)/sqrt(nrow(raw)), n_reps, length(coverages)))
cat(sprintf("alpha = %.2f, full effect VE = %.2f at f = 1\n", alpha_vac, 1 - alpha_vac))

dir.create("output", showWarnings = FALSE)
fwrite(res, "output/rct_cir_vs_ve.csv")

long <- melt(res, id.vars = "coverage", measure.vars = c("VE_trial", "VE_flip"),
             variable.name = "estimand", value.name = "VE")
long[, estimand := factor(estimand, levels = c("VE_trial", "VE_flip"),
                          labels = c("Trial contrast (1 - CIR)",
                                     "Causal flip-VE"))]
p <- ggplot(long, aes(coverage, VE, colour = estimand)) +
  geom_hline(yintercept = 1 - alpha_vac, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_x_continuous(labels = scales::percent) +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(x = "vaccine coverage", y = "VE",
       title = if (thresh > 0)
         "Threshold vaccine effect: trial contrast vs causal VE across coverage"
       else
         "Contact-dependent effect: trial contrast vs causal VE across coverage",
       subtitle = sprintf("N = %d, dashed = 1 - alpha = %.2f", N, 1 - alpha_vac))
ggsave("output/rct_cir_vs_ve.png", p, width = 8, height = 5, dpi = 130)
cat("Wrote output/rct_cir_vs_ve.csv and .png\n")
