# Two local mechanisms that make the observed CIR differ from the causal VE
# ---------------------------------------------------------------------------
# Earlier sweeps found CIR = 1 - AR_vac/AR_unvac to be an unbiased estimate of
# the causal own-outcome flip-VE at every coverage, under a fixed effect, a
# smooth contact-dependent effect and a threshold effect, on directed and
# symmetric graphs (pooled gap -0.004 +/- 0.007). Randomisation makes the trial
# identify the direct effect; interference in the vaccine's own MECHANISM does
# not break that, because the spillover from flipping i lands on other people's
# outcomes, which are counted in their own flips.
#
# Both mechanisms here break that, locally and without fine-tuning:
#
#   A  trans_tau < 1  -- coverage-dependent TRANSMISSIBILITY. A vaccinated
#      person in a well-vaccinated neighbourhood is also less infectious. Now
#      flipping i directly lowers the risk i poses to its contacts, so the
#      spillover lands INSIDE i's own flip and the flip-VE sees it while the
#      vaccinated-vs-control contrast does not.
#
#   B  spill > 0      -- NEIGHBOURHOOD IMMUNITY. Unvaccinated people are also
#      protected by living among vaccinated contacts, at `spill` of the full
#      effect. This contaminates the control arm: a control's outcome depends
#      on other people's treatment, so the CIR is biased toward 1 (SUTVA
#      violated) while the flip-VE is not.
#
# Both keep the effect LOCAL (it depends on a node's own neighbourhood) and
# SMOOTH (linear in local coverage) -- no knife-edge on a population total.
#
# beta is calibrated per configuration so the unvaccinated attack rate is
# ~target_ar throughout, so nothing here is a saturation artefact.
#
#   Rscript cir_ve_spillover.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()

N         <- 1000
mean_k    <- 6
alpha_vac <- 0.1          # full effect VE = 0.9 at f = 1
target_ar <- 0.5
gamma     <- 1
t_star    <- 15
init_I    <- 5
n_rep     <- 2500
n_rep_ve  <- 400
n_flip    <- 120
n_reps    <- 3            # independent trials per cell
cores     <- 8
coverages <- c(0.25, 0.45, 0.65, 0.85)

arms <- list(
  list(id = "baseline",        spill = 0.0, tau = 1.0),
  list(id = "A transmissible", spill = 0.0, tau = 0.2),
  list(id = "B spillover",     spill = 0.6, tau = 1.0),
  list(id = "A + B",           spill = 0.6, tau = 0.2))

set.seed(1)
adj <- sample_pareto_adj(N, alpha = 3, mean_k = mean_k, seed = 1)
csr <- adj_to_csr(adj = adj)
tp  <- seq(1, t_star, 1)
cat(sprintf("N=%d mean_deg=%.2f alpha=%.2f  target AR=%.2f  %d reps/cell\n\n",
            N, mean(adj$degree), alpha_vac, target_ar, n_reps))

sim <- function(vac, sus, b, nsim, seed) {
  run_stoch_network_events(
    beta = b, N = N, susceptibility = sus, t = t_star, vac = vac, csr = csr,
    gamma = gamma, timepoints = tp, I_ini = init_I, n_sim = nsim, seed = seed,
    k_mean = mean_k, cores = cores, return_times = TRUE,
    transmissibility = attr(sus, "trans"))
}

cell <- function(cov, arm, rep) {
  set.seed(3000 + round(1000 * cov) + 7919L * rep)
  vac   <- sample.int(N, round(cov * N))
  unvac <- setdiff(seq_len(N), vac)
  mk <- function(v) vacfrac_susceptibility(v, alpha_vac, adj = adj,
                                           vac_frac_ref = 1, vac_frac_power = 1,
                                           spill = arm$spill, trans_tau = arm$tau)
  sus <- mk(vac)

  lo <- log(0.2); hi <- log(120)
  for (it in 1:15) {
    m <- (lo + hi) / 2
    ar <- mean(colMeans(sim(vac, sus, exp(m), 400, 3L + rep) <= t_star)[unvac])
    if (ar < target_ar) lo <- m else hi <- m
  }
  b <- exp((lo + hi) / 2)

  P    <- colMeans(sim(vac, sus, b, n_rep, 7L + 100L * rep) <= t_star)
  ar_v <- mean(P[vac]); ar_u <- mean(P[unvac])

  # Causal flip-VE by re-simulation: flip i, take i's own outcome. CRN across
  # the factual and flipped runs so the contrast is not swamped by noise.
  set.seed(11 + rep)
  fl_v <- sample(vac,   min(round(n_flip * cov), length(vac)))
  fl_u <- sample(unvac, min(n_flip - length(fl_v), length(unvac)))
  sd_i <- 500L + rep
  p_i  <- function(i, v2) colMeans(sim(v2, mk(v2), b, n_rep_ve, sd_i) <= t_star)[i]
  base <- colMeans(sim(vac, sus, b, n_rep_ve, sd_i) <= t_star)
  cf_u <- if (length(fl_u)) vapply(fl_u, function(i) p_i(i, sort(c(vac, i))), 0) else numeric(0)
  cf_v <- if (length(fl_v)) vapply(fl_v, function(i) p_i(i, setdiff(vac, i)),  0) else numeric(0)

  num   <- sum(base[vac])   + length(unvac) * mean(cf_u)
  denom <- sum(base[unvac]) + length(vac)   * mean(cf_v)

  data.table(arm = arm$id, coverage = cov, rep = rep, beta = b,
             AR_unvac = ar_u, AR_vac = ar_v,
             VE_trial = 1 - ar_v / ar_u, VE_flip = 1 - num / denom)
}

raw <- rbindlist(lapply(arms, function(a)
  rbindlist(lapply(coverages, function(cv)
    rbindlist(lapply(seq_len(n_reps), function(rp) cell(cv, a, rp)))))))
raw[, gap := VE_flip - VE_trial]
fwrite(raw, "output/cir_ve_spillover_reps.csv")

res <- raw[, .(beta = mean(beta), AR_unvac = mean(AR_unvac),
               VE_trial = mean(VE_trial), VE_flip = mean(VE_flip),
               gap = mean(gap), se = sd(gap) / sqrt(.N)),
           by = .(arm, coverage)]
res[, t := gap / se]

cat("=== trial contrast vs causal flip-VE ===\n")
print(res[, .(arm, coverage, beta = round(beta, 2), AR_unvac = round(AR_unvac, 3),
              VE_trial = round(VE_trial, 3), VE_flip = round(VE_flip, 3),
              gap = round(gap, 3), se = round(se, 3), t = round(t, 1))])
cat("\npooled gap by arm (|t| > 2 means the CIR is measurably NOT the causal VE):\n")
print(raw[, .(gap = round(mean(gap), 4), se = round(sd(gap) / sqrt(.N), 4),
              t = round(mean(gap) / (sd(gap) / sqrt(.N)), 1)), by = arm])
fwrite(res, "output/cir_ve_spillover.csv")

long <- melt(res, id.vars = c("arm", "coverage"),
             measure.vars = c("VE_trial", "VE_flip"),
             variable.name = "estimand", value.name = "VE")
long[, estimand := factor(estimand, levels = c("VE_trial", "VE_flip"),
                          labels = c("Trial contrast (1 - CIR)", "Causal flip-VE"))]
p <- ggplot(long, aes(coverage, VE, colour = estimand)) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  facet_wrap(~ arm) +
  scale_x_continuous(labels = scales::percent) +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = "vaccine coverage", y = "VE",
       title = "When does the observed CIR stop measuring the causal VE?",
       subtitle = "attack rate held at 0.5 throughout; A = coverage-dependent transmissibility, B = neighbourhood immunity")
ggsave("output/cir_ve_spillover.png", p, width = 9, height = 6, dpi = 130)
cat("\nWrote output/cir_ve_spillover.csv and .png\n")
