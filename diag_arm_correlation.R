# Why the trial contrast is precise in some models and useless in others
# ---------------------------------------------------------------------------
# VE is a RATIO, so what matters for its precision is not how variable the
# epidemic is but how much of that variability is SHARED between the two arms.
# Variation common to both cancels; variation independent between them does not.
#
# The four models here make that a controlled experiment, because two of them
# have almost identical level variability and opposite arm correlation:
#
#                    cv(C1)   cor(C1,C2)   sd(1 - CIR)
#   linear            0.061      0.008         0.059
#   SIR (mixed)       0.471      0.971         0.148
#   network pa=2      0.630      0.985         0.255
#   segregated        0.627      0.001        22.3
#
# The network and the segregated model differ by 0.3% in cv(C1) and by a factor
# of ~90 in the precision of the estimand. The only difference between them is
# whether the arms share an epidemic.
#
# Mechanically: in the segregated (two-cluster) design each arm is its own
# independent epidemic that can fizzle or take off. When the control cluster
# fizzles and the vaccinated one does not, C2/C1 explodes and 1 - CIR goes
# hugely negative. 5% of replicates have C1 = 0 outright, where the ratio does
# not exist at all.
#
# This is the cluster-RCT design effect in its sharpest form: not a factor of
# 1 + (m-1)*rho, but a ratio estimand with essentially no precision despite each
# arm being perfectly well measured on its own.
#
# It is also a calibration check on the kernel posterior. A predictive spread
# this wide should give a flat likelihood and hence a very wide posterior for
# alpha and VE under sir_segregated. If that model's reported interval comes
# back comparable to the others, the posterior is not capturing it.
#
# Each model is at its own calibration to the shared target (C1 ~ 160, C2 ~ 80
# at t* = 8), so the comparison is at matched data rather than matched
# parameters.
#
#   Rscript diag_arm_correlation.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table); library(ggplot2)})
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()
`%||%` <- function(a, b) if (is.null(a)) b else a

# materialise_cfg / build_simulator are defined in run_fit_array.R; pull just
# those two in rather than sourcing it (which would launch the whole array).
src <- readLines("run_fit_array.R")
for (fn in c("materialise_cfg", "build_simulator")) {
  a <- grep(paste0("^", fn, " <- function"), src)
  b <- a + which(src[a:length(src)] == "}")[1] - 1
  eval(parse(text = paste(src[a:b], collapse = "\n")))
}

N <- 800; t_star <- 8; gamma <- 1; dt <- 0.05; B <- 4000
set.seed(11); vac0 <- sample.int(N, 400); unvac0 <- setdiff(seq_len(N), vac0)

seg <- {
  cfg <- materialise_cfg(list(
    model_type = "sir_segregated", sim_type = "sir_multisite",
    N_cont = 400, N_vac = 400, t_star = t_star, gamma = gamma,
    I_ini_2g = c(3, 3), dt = dt, n_sites = 2L, site_icc = 1,
    allocation_seed = 1L, inner_cores = 8))
  o <- build_simulator(cfg)(1.630, 0.8358, B, seed = 2024); setDT(o)
  data.table(model = "segregated (2-cluster)", C1 = o$C1, C2 = o$C2)
}
lin <- {
  r <- run_stoch_linear_dust(beta = 0.0637, N = c(400, 400),
        susceptibility = c(1, 0.4392), t = t_star, dt = dt,
        timepoints = t_star, n_sim = B, cores = 8)
  setDT(r); f <- r[time == t_star]
  data.table(model = "linear", C1 = f$C1, C2 = f$C2)
}
sir <- {
  r <- run_stoch_cd_dust(matrix(1, 2, 2), beta = 2.1045, N = c(400, 400),
        t = t_star, I_ini = c(3, 3), susceptibility = c(1, 0.4147),
        gamma = gamma, dt = dt, timepoints = t_star, n_sim = B, cores = 8,
        seed = 2024)
  setDT(r); f <- r[time == t_star]
  data.table(model = "SIR (mixed)", C1 = f$C1, C2 = f$C2)
}
net <- {
  set.seed(1); cm <- get_conact_matrix_pl(N, alpha = 2, mean_k = 6)
  csr <- adj_to_csr(contact_matrix = cm, adj = contact_matrix_to_adj(cm))
  sus <- rep(1, N); sus[vac0] <- 0.2820
  inf <- run_stoch_network_events(beta = 2.723, N = N, susceptibility = sus,
          t = t_star, vac = vac0, csr = csr, gamma = gamma,
          timepoints = t_star, I_ini = 3, n_sim = B, seed = 2024,
          k_mean = 6, cores = 8, return_times = TRUE)
  h <- inf <= t_star
  data.table(model = "network pa=2", C1 = rowSums(h[, unvac0, drop = FALSE]),
             C2 = rowSums(h[, vac0, drop = FALSE]))
}

d <- rbindlist(list(lin, sir, net, seg))
d[, `:=`(tot = C1 + C2, ve_raw = fifelse(C1 > 0, 1 - C2 / C1, NA_real_))]
fwrite(d, "output/arm_correlation.csv")

out <- d[, {
  ok <- !is.na(ve_raw)
  .(mean_C1 = round(mean(C1), 1), mean_C2 = round(mean(C2), 1),
    cv_C1 = round(sd(C1) / mean(C1), 3),
    cor_arms = round(cor(C1, C2), 4),
    p_C1_zero = round(mean(C1 == 0), 3),
    ratio_mean = round(mean(ve_raw[ok]), 3),
    sd_ratio = round(sd(ve_raw[ok]), 4),
    q10 = round(quantile(ve_raw[ok], .1), 3),
    q90 = round(quantile(ve_raw[ok], .9), 3))
}, by = model]
cat(sprintf("%d replicate trials per model, each at its own calibration to C1 ~ 160 / C2 ~ 80\n\n", B))
print(out)
cat("\nsd_ratio is the spread of the trial's OWN VE estimate (1 - CIR) across\n")
cat("replicate trials. cv_C1 is how variable the epidemic size is.\n")
cat("\nnetwork vs segregated: cv_C1 differs by 0.3%, sd_ratio by ~90x.\n")
cat("The only difference is whether the two arms share an epidemic.\n")

# Scatter: the arms move together, or they do not.
p <- ggplot(d, aes(C1, C2)) +
  geom_point(alpha = 0.12, size = 0.6) +
  facet_wrap(~ model, scales = "free") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = "C1 (unvaccinated cases)", y = "C2 (vaccinated cases)",
       title = "Whether the arms share an epidemic decides the estimand's precision",
       subtitle = "a tight diagonal means the epidemic-size variation cancels in the ratio")
ggsave("output/arm_correlation.png", p, width = 9, height = 7, dpi = 140)
cat("\nWrote output/arm_correlation.{csv,png}\n")
