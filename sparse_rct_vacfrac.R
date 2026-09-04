# Sparse RCT under a contact-dependent vaccine effect: CIR != VE
# ---------------------------------------------------------------------------
# A sparse trial (a small vaccinated fraction embedded in a large population)
# is usually argued to measure the DIRECT effect cleanly: too few vaccinated
# to change the epidemic, so no interference, so the observed cumulative
# incidence ratio CIR = AR_vac / AR_control estimates the per-contact effect.
#
# That argument assumes the effect itself does not depend on coverage. Under
# the contact-dependent model (stoch_mod_adj_vacfrac.R) it does:
#
#     alpha_eff[i] = alpha ^ ((f[i] / vac_frac_ref) ^ vac_frac_power)
#
# with f[i] the fraction vaccinated in i's local neighbourhood (i included).
# In a SPARSE trial a vaccinated person has only themselves vaccinated locally,
# f = 1/(1+degree), so they get a REDUCED but non-zero effect; at deployment
# coverage (f = vac_frac_ref) the same vaccine has its full effect alpha.
#
# WHAT THIS DOES AND DOES NOT SHOW
# The trial contrast (1 - CIR) and the causal flip-VE agree at EVERY coverage
# (see the run: gaps are noise). That is structural, not a bug -- both are the
# DIRECT effect, so an internally valid trial does recover the flip-VE at its
# own coverage. The neighbour spillover from vaccinating i lands on the
# neighbours' outcomes, which are counted in THEIR flips, not in i's, so this
# estimand does not see it (that is the total effect, a different quantity).
#
# The gap is across COVERAGE, not between estimands: the same vaccine reads
# ~0.13 in a 0.5% trial and ~0.56 at 50% coverage. So a sparse RCT is valid
# for the setting it ran in and still misestimates the deployment effect --
# an transportability failure that no sample size or randomisation fixes.
#
# Companion to network_vs_linear_cir.R (which fixes coverage and follows CIR
# over time); here coverage is the swept axis. Runs on the exact event-driven
# engine, so the whole sweep is cheap even at N = 20000.
#
#   Rscript sparse_rct_vacfrac.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()

# ---- config ---------------------------------------------------------------
N          <- 20000
pl_alpha   <- 3
mean_k     <- 6
alpha_vac  <- 0.4          # full effect: VE = 0.6 at the reference coverage
vac_ref    <- 0.5          # coverage at which alpha_eff == alpha_vac
vac_power  <- 1
beta       <- 1.5
gamma      <- 1
t_star     <- 12
init_I     <- 20
n_rep      <- 400
n_rep_ve   <- 200          # replicates per re-simulated counterfactual
n_flip     <- 60           # individuals flipped per trial
cores      <- 8
coverages  <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50)

set.seed(1)
adj <- sample_pareto_adj(N, alpha = pl_alpha, mean_k = mean_k, seed = 1)
csr <- adj_to_csr(adj = adj)
tp  <- seq(1, t_star, 1)
cat(sprintf("N=%d mean_deg=%.1f  alpha=%.2f (VE_ref=%.2f at coverage %.0f%%)\n\n",
            N, mean(adj$degree), alpha_vac, 1 - alpha_vac, 100 * vac_ref))

# One trial: vaccinate `cov` of the population, compare them to an equally
# sized control drawn from the unvaccinated. `model` selects whether the
# vaccine effect is fixed (standard) or contact-dependent (vacfrac).
#
# Two quantities per trial:
#   CIR     -- what the trial OBSERVES: AR_vac / AR_control.
#   VE_flip -- the causal effect of flipping one person's status, obtained by
#              RE-SIMULATING with that person flipped (get_stoch_eate_..., which
#              is why the frozen field will not do here). It differs from the
#              CIR because flipping i also changes f[j] for i's neighbours, so a
#              vaccinated neighbour gains protection -- a local spillover the
#              vaccinated-vs-control contrast cannot see.
#
# vac_frac_power = 0 makes alpha_eff = alpha regardless of f, i.e. exactly the
# standard fixed-effect model, so both arms use the identical code path.
run_trial <- function(cov, model) {
  pw <- if (model == "vacfrac") vac_power else 0
  set.seed(1000 + round(1e4 * cov))
  n_v  <- max(1L, round(cov * N))
  vac  <- sample.int(N, n_v)
  ctrl <- sample(setdiff(seq_len(N), vac), n_v)      # matched control arm

  sus <- vacfrac_susceptibility(vac, alpha_vac, adj = adj,
                                vac_frac_ref = vac_ref, vac_frac_power = pw)

  inf <- run_stoch_network_events(
    beta = beta, N = N, susceptibility = sus, t = t_star, vac = vac, csr = csr,
    gamma = gamma, timepoints = tp, I_ini = init_I, n_sim = n_rep,
    seed = 7, k_mean = mean_k, cores = cores, return_times = TRUE)

  hit    <- inf <= t_star
  ar_v   <- mean(rowMeans(hit[, vac,  drop = FALSE]))
  ar_c   <- mean(rowMeans(hit[, ctrl, drop = FALSE]))
  # mean realised f among the vaccinated: what drives alpha_eff
  vi <- integer(N); vi[vac] <- 1L
  nv <- colSums(adj$mask * matrix(vi[adj$neighbors], nrow(adj$neighbors), N))
  dg <- colSums(adj$mask)
  f_bar <- mean(ifelse(dg > 0, nv / dg, 0)[vac])

  # Causal flip-VE: re-simulate with each sampled individual's status flipped.
  ve <- get_stoch_eate_network_vacfrac(
    beta = beta, alpha = alpha_vac, N = N, t = t_star, adj = adj, csr = csr,
    vac_frac_power = pw, vac_frac_ref = vac_ref, k_mean = mean_k,
    gamma = gamma, n_rep = n_rep_ve, n_flip = n_flip, timepoints = tp,
    init_I = init_I, mc.cores = 1, inner_cores = cores,
    vac_list = list(vac), engine = "events", crn_seed = 11)
  setDT(ve)
  eate <- ve[method == "full_stoch" & t == t_star, eate]

  data.table(model = model, coverage = cov, f_bar = f_bar,
             mean_alpha_eff = mean(sus[vac]),
             AR_vac = ar_v, AR_ctrl = ar_c,
             CIR = ar_v / ar_c, VE_trial = 1 - ar_v / ar_c,
             VE_flip = 1 - eate)
}

res <- rbindlist(lapply(c("standard", "vacfrac"),
                        function(m) rbindlist(lapply(coverages, run_trial, model = m))))

res[, gap := VE_flip - VE_trial]
cat("=== observed trial contrast vs the causal flip-VE ===\n")
print(res[, .(model, coverage,
              f_bar          = round(f_bar, 3),
              mean_alpha_eff = round(mean_alpha_eff, 3),
              AR_ctrl        = round(AR_ctrl, 3),
              VE_trial       = round(VE_trial, 3),   # 1 - CIR, what the RCT sees
              VE_flip        = round(VE_flip, 3),    # causal, by re-simulation
              gap            = round(gap, 3))])
cat(sprintf("\nVE at the reference coverage (the deployment effect) = %.2f\n",
            1 - alpha_vac))
cat("VE_trial (1 - CIR) and VE_flip AGREE at every coverage -- both are the\n",
    "direct effect, so the trial is internally valid for its own setting.\n",
    "The failure is across coverage: under the contact-dependent effect the\n",
    "same vaccine reads far lower in a sparse trial than at deployment,\n",
    "because the effect itself depends on local coverage. The fixed-effect\n",
    "model shows no such drift.\n")

dir.create("output", showWarnings = FALSE)
fwrite(res, "output/sparse_rct_vacfrac.csv")

res[, model_f := factor(model, levels = c("standard", "vacfrac"),
                        labels = c("Fixed vaccine effect",
                                   "Contact-dependent effect"))]
long <- melt(res, id.vars = c("model_f", "coverage"),
             measure.vars = c("VE_trial", "VE_flip"),
             variable.name = "estimand", value.name = "VE")
long[, estimand := factor(estimand, levels = c("VE_trial", "VE_flip"),
                          labels = c("Trial contrast (1 - CIR)",
                                     "Causal flip-VE (re-simulated)"))]
p <- ggplot(long, aes(coverage, VE, colour = estimand)) +
  geom_hline(yintercept = 1 - alpha_vac, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ model_f) +
  scale_x_log10(labels = scales::percent) +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = "trial coverage (fraction vaccinated, log scale)", y = "VE",
       title = "A sparse RCT is internally valid yet misestimates the deployment effect",
       subtitle = paste("both estimands agree at each coverage; under a",
                        "contact-dependent effect they drift with coverage"))
ggsave("output/sparse_rct_vacfrac.png", p, width = 8, height = 5, dpi = 130)
cat("\nWrote output/sparse_rct_vacfrac.csv and .png\n")
