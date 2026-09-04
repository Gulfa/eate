# The observed CIR and the causal VE can differ by an ARBITRARY amount
# ---------------------------------------------------------------------------
# Everything else in this repo asks how big the gap between the trial contrast
#     CIR = AR_vac / AR_unvac
# and the causal estimand
#     EATE = E_i[Y_i(V_i = 1)] / E_i[Y_i(V_i = 0)]
# gets under some plausible mechanism. Answer, empirically: small under a
# homogeneous mechanism at any coverage (~0.01), a few hundredths under
# degree heterogeneity, ~0.1 on a heavy-tailed network with a
# contact-dependent effect.
#
# This script makes the complementary point: those are facts about the
# mechanisms we chose, not about the estimands. With NO constraint on what a
# vaccine may do, the gap is unbounded. Take
#
#     alpha_i = alpha1  if the number vaccinated is ODD
#               alpha2  if it is EVEN
#
# Flipping any single person changes that count by one, so it flips the parity,
# so it changes EVERYONE's susceptibility. Every counterfactual Y_i(V_i = 1 - v_i)
# is therefore evaluated in a different epidemic from the factual one, and the
# gap is essentially VE(alpha1) - VE(alpha2): choose those two far apart and the
# gap is as large as you like, in either direction.
#
# Nothing here is subtle -- the point is precisely that it is not. The mechanism
# is absurd as biology, and no experiment can rule it out from outcome data,
# because the factual world only ever shows one parity. So bounding the gap
# between what a trial measures and what we want to know is a MODELLING
# commitment ("effects are local, smooth, individually mediated"), not something
# the data or the randomisation can deliver.
#
#   Rscript parity_ve_unbounded.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R")

# ---- config ---------------------------------------------------------------
N        <- 400            # total; half vaccinated
f        <- 0.5
beta     <- 2
gamma    <- 1
t_star   <- 12
I_ini    <- c(5, 5)
dt       <- 0.05
n_rep    <- 4000
alpha1   <- 0.2            # susceptibility when the vaccinated count is ODD
alpha2s  <- c(0.2, 0.4, 0.6, 0.8, 1.0)   # ... and when it is EVEN
mm       <- matrix(1, 2, 2)

# Attack rates in the two arms for a given alpha, everything else fixed.
arm_ar <- function(alpha, n_vac, seed) {
  n_unvac <- N - n_vac
  raw <- run_stoch_cd_dust(mm, beta = beta, N = c(n_unvac, n_vac), t = t_star,
                           I_ini = I_ini, susceptibility = c(1, alpha),
                           gamma = gamma, dt = dt, timepoints = t_star,
                           n_sim = n_rep, cores = 8, seed = seed)
  setDT(raw); fin <- raw[time == t_star]
  c(unvac = mean(fin$C1) / n_unvac, vac = mean(fin$C2) / n_vac)
}

# Put the FACTUAL count on the odd side, so the factual world always runs at
# alpha1 and is identical in every row: the trial sees exactly the same data
# whatever alpha2 is. Only the counterfactuals move.
n_vac_fac <- round(f * N); if (n_vac_fac %% 2 == 0) n_vac_fac <- n_vac_fac + 1L

res <- rbindlist(lapply(alpha2s, function(alpha2) {
  a_of <- function(n) if (n %% 2 == 1) alpha1 else alpha2

  # FACTUAL world: n_vac_fac vaccinated, so alpha = a_of(n_vac_fac).
  fac <- arm_ar(a_of(n_vac_fac), n_vac_fac, seed = 11)
  cir <- fac[["vac"]] / fac[["unvac"]]

  # COUNTERFACTUALS. Flipping any one person changes the count by 1, flipping
  # the parity, so every counterfactual world runs at the OTHER alpha.
  #   - an unvaccinated i, if vaccinated: count + 1
  #   - a vaccinated i, if unvaccinated:  count - 1
  cf_up <- arm_ar(a_of(n_vac_fac + 1), n_vac_fac + 1, seed = 12)
  cf_dn <- arm_ar(a_of(n_vac_fac - 1), n_vac_fac - 1, seed = 13)

  # EATE = E_i[Y_i(1)] / E_i[Y_i(0)], each individual contributing to both:
  # observed where their status matches the arm, counterfactual where it does not.
  n_v <- n_vac_fac; n_u <- N - n_vac_fac
  num <- n_v * fac[["vac"]]   + n_u * cf_up[["vac"]]
  den <- n_u * fac[["unvac"]] + n_v * cf_dn[["unvac"]]

  data.table(alpha1 = alpha1, alpha2 = alpha2,
             AR_unvac = fac[["unvac"]], AR_vac = fac[["vac"]],
             VE_trial = 1 - cir, VE_eate = 1 - num / den)
}))
res[, gap := VE_eate - VE_trial]

cat(sprintf("N=%d, %.0f%% vaccinated, alpha = %.2f when the vaccinated count is ODD\n\n",
            N, 100 * f, alpha1))
print(res[, .(alpha2, AR_unvac = round(AR_unvac, 3), AR_vac = round(AR_vac, 3),
              VE_trial = round(VE_trial, 3), VE_eate = round(VE_eate, 3),
              gap = round(gap, 3))])
cat("\nThe trial contrast is unchanged across rows -- the factual world is the same.\n",
    "Only the counterfactual parity differs, and the gap tracks alpha2 with no\n",
    "bound: the observed data cannot see it, because it never shows both parities.\n")

dir.create("output", showWarnings = FALSE)
fwrite(res, "output/parity_ve_unbounded.csv")

long <- melt(res, id.vars = "alpha2", measure.vars = c("VE_trial", "VE_eate"),
             variable.name = "estimand", value.name = "VE")
long[, estimand := factor(estimand, levels = c("VE_trial", "VE_eate"),
                          labels = c("Trial contrast (1 - CIR)", "Causal EATE"))]
p <- ggplot(long, aes(alpha2, VE, colour = estimand)) +
  geom_line(linewidth = 1) + geom_point(size = 2.4) +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(x = expression(alpha[2]~"(susceptibility when the vaccinated count is even)"),
       y = "VE",
       title = "Without modelling constraints the trial contrast and the causal VE are unrelated",
       subtitle = sprintf("alpha = %.2f when odd; the factual world -- and so the trial -- is identical in every case",
                          alpha1))
ggsave("output/parity_ve_unbounded.png", p, width = 8, height = 5, dpi = 130)
cat("Wrote output/parity_ve_unbounded.csv and .png\n")
