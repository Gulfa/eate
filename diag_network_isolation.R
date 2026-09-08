# How often does sample_pareto_adj produce an unusable contact network?
# ---------------------------------------------------------------------------
# At heavy tails the Pareto propensities concentrate in a few hubs, so a large
# share of nodes end up with degree 0 and can never be infected. That is what
# made some network configs unfittable: seed 10 at pl_alpha = 1.2 had 27%
# isolated nodes and a mean C1 of 0.8 against a target of 160, so the kernel
# likelihood was 0 everywhere and the fit pinned at the parameter bound.
#
# This sweeps pl_alpha from 1 to 3 at mean_k = 6 and 12, with 20 network
# realisations each, and reports the isolated fraction (and the realised mean
# degree, which falls well short of mean_k at heavy tails).
#
#   Rscript diag_network_isolation.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R")

N        <- 800
n_seed   <- 20
pl_alphas <- c(1, 1.2, 1.5, 2, 2.5, 3)
mean_ks  <- c(6, 12)

res <- rbindlist(lapply(mean_ks, function(mk)
  rbindlist(lapply(pl_alphas, function(pa)
    rbindlist(lapply(seq_len(n_seed), function(sd_) {
      adj <- sample_pareto_adj(N, alpha = pa, mean_k = mk, seed = sd_)
      d   <- adj$degree
      data.table(mean_k = mk, pl_alpha = pa, seed = sd_,
                 iso = mean(d == 0), deg = mean(d), maxdeg = max(d))
    }))))))

fwrite(res, "output/network_isolation.csv")

summ <- res[, .(iso_med = median(iso), iso_lo = min(iso), iso_hi = max(iso),
                deg_med = median(deg), maxdeg_med = median(maxdeg),
                n_bad = sum(iso > 0.15)),
            by = .(mean_k, pl_alpha)][order(mean_k, pl_alpha)]
cat("=== isolated fraction over", n_seed, "realisations (N =", N, ") ===\n")
print(summ[, .(mean_k, pl_alpha,
               iso_med = round(iso_med, 3), iso_min = round(iso_lo, 3),
               iso_max = round(iso_hi, 3), deg_med = round(deg_med, 2),
               maxdeg = maxdeg_med, n_over_15pct = n_bad)])

res[, mk := factor(mean_k, levels = mean_ks,
                   labels = paste0("mean_k = ", mean_ks))]
p <- ggplot(res, aes(factor(pl_alpha), iso, colour = mk)) +
  geom_hline(yintercept = 0.15, linetype = "dashed", colour = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.55, position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.7),
             alpha = 0.55, size = 1.5) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(x = "Pareto exponent (pl_alpha)", y = "isolated nodes (degree 0)",
       title = sprintf("Heavy-tailed contact networks strand a large share of nodes (N = %d, %d realisations)",
                       N, n_seed),
       subtitle = "dashed line = 15%, above which configs in this repo failed to fit")
ggsave("output/network_isolation.png", p, width = 9, height = 5.5, dpi = 140)

# Realised mean degree vs the requested mean_k -- the other half of the story.
p2 <- ggplot(res, aes(factor(pl_alpha), deg, colour = mk)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.7),
             alpha = 0.55, size = 1.5) +
  geom_hline(yintercept = mean_ks, linetype = "dotted", colour = "grey55") +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(x = "Pareto exponent (pl_alpha)", y = "realised mean degree",
       title = "Realised mean degree falls short of the requested mean_k at heavy tails",
       subtitle = "dotted lines = requested mean_k (6 and 12)")
ggsave("output/network_isolation_degree.png", p2, width = 9, height = 5.5, dpi = 140)

cat("\nWrote output/network_isolation.{csv,png} and network_isolation_degree.png\n")
