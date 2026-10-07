# What a trial measures vs what the vaccine actually does
# ---------------------------------------------------------------------------
# Two axes:
#   x  1 - CIR, the VE the trial would report from its own arm contrast
#   y  infections averted per 1000, against the same population unvaccinated
# and the curve traced out as vaccine coverage is swept.
#
# Three mechanisms, same beta and same per-contact effect size, so the only
# thing differing is HOW the vaccine acts:
#
#   linear          no transmission feedback at all. Vaccination cannot
#                   protect anyone but the recipient.
#   SIR, blocks     reduces SUSCEPTIBILITY. Individual protection plus herd
#   infection       immunity, so the trial sees something and the population
#                   gets more than the trial suggests.
#   SIR, blocks     reduces ONWARD TRANSMISSION only, leaving susceptibility
#   transmission    untouched. Vaccinated and unvaccinated mix and so face the
#                   same force of infection: the arm contrast is ~1 and the
#                   trial reports VE ~ 0, while the population effect is large.
#
# That last one is the point of the figure. A trial contrast is a statement
# about one individual's risk relative to another's in the SAME epidemic; it is
# blind to an effect that acts entirely through the epidemic itself.
#
# The no-vaccination baseline is the same model with the vaccine parameter set
# to 1, rather than a zero-size arm, so the group structure and seeding are
# identical across the comparison.
#
#   Rscript diag_indirect_effects.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R")

N <- 2000; t_star <- 12; gamma <- 1; dt <- 0.05
beta_sir <- 2.0           # R0 = 2 unvaccinated
beta_lin <- 0.075         # matched so the unvaccinated attack rate is similar
I_tot    <- 10
eff      <- 0.4           # susceptibility multiplier, or transmissibility one
n_sim    <- 3000
covs     <- c(0.05, seq(0.1, 0.9, by = 0.1), 0.95)

# One run of the two-group model at coverage f. `mode` picks which channel the
# vaccine acts through; v = 1 is the no-vaccination baseline.
run_at <- function(f, mode, v, seed = 11) {
  n_v <- max(1L, round(f * N)); n_u <- N - n_v
  Ngr <- c(n_u, n_v)
  Ii  <- .spread_seeds(I_tot, Ngr)
  if (mode == "linear") {
    r <- run_stoch_linear_dust(beta = beta_lin, N = Ngr,
          susceptibility = c(1, v), t = t_star, dt = dt,
          timepoints = t_star, n_sim = n_sim, cores = 8)
  } else {
    sus <- if (mode == "susceptibility") c(1, v) else c(1, 1)
    trn <- if (mode == "transmission")   c(1, v) else c(1, 1)
    r <- run_stoch_cd_dust(matrix(1, 2, 2), beta = beta_sir, N = Ngr,
          t = t_star, I_ini = Ii, susceptibility = sus,
          transmissibility = trn, gamma = gamma, dt = dt,
          timepoints = t_star, n_sim = n_sim, cores = 8, seed = seed)
  }
  setDT(r); fin <- r[time == t_star]
  list(AR_u = mean(fin$C1) / n_u, AR_v = mean(fin$C2) / n_v,
       total = mean(fin$C1 + fin$C2))
}

modes <- c(linear = "linear", susceptibility = "susceptibility",
           transmission = "transmission")
res <- rbindlist(lapply(names(modes), function(m) {
  base <- run_at(0.5, m, v = 1)$total          # same model, vaccine switched off
  rbindlist(lapply(covs, function(f) {
    o <- run_at(f, m, v = eff)
    data.table(mode = m, coverage = f,
               AR_unvac = o$AR_u, AR_vac = o$AR_v,
               CIR = o$AR_v / o$AR_u,
               ve_trial = 1 - o$AR_v / o$AR_u,
               total = o$total, baseline = base,
               averted_per1k = (base - o$total) / N * 1000)
  }))
}))
fwrite(res, "output/indirect_effects.csv")

lab <- c(linear = "Linear (no transmission feedback)",
         susceptibility = "SIR, blocks infection",
         transmission = "SIR, blocks onward transmission only")
res[, model := factor(lab[mode], levels = unname(lab))]

cat(sprintf("N = %d, t* = %d, vaccine parameter = %.2f, baseline = same model unvaccinated\n\n",
            N, t_star, eff))
print(res[, .(model = mode, coverage, ve_trial = round(ve_trial, 3),
              averted_per1k = round(averted_per1k, 1))])

p <- ggplot(res, aes(ve_trial, averted_per1k, colour = model)) +
  geom_path(linewidth = 0.9) +
  geom_point(size = 1.9) +
  geom_text(data = res[coverage %in% c(0.1, 0.5, 0.9)],
            aes(label = sprintf("%.0f%%", 100 * coverage)),
            hjust = -0.35, vjust = 0.4, size = 3, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(x = "VE the trial would report  (1 - CIR)",
       y = "infections averted per 1000, vs no vaccination",
       title = "What the trial measures against what the vaccine does",
       subtitle = "each curve is one mechanism, swept over coverage (labelled 10 / 50 / 90%)")
ggsave("output/indirect_effects.png", p, width = 9, height = 6, dpi = 140)

p2 <- ggplot(melt(res, id.vars = c("model", "coverage"),
                  measure.vars = c("ve_trial", "averted_per1k")),
             aes(coverage, value, colour = model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  facet_wrap(~ variable, scales = "free_y",
             labeller = as_labeller(c(ve_trial = "VE the trial reports (1 - CIR)",
                                      averted_per1k = "infections averted per 1000"))) +
  scale_x_continuous(labels = scales::percent) +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = "vaccine coverage", y = NULL,
       title = "The same two quantities against coverage")
ggsave("output/indirect_effects_by_coverage.png", p2, width = 10, height = 5, dpi = 140)
cat("\nWrote output/indirect_effects.csv and the two figures\n")
