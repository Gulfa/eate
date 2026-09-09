# How much of a wide beta posterior actually reaches the AVE?
#
# SCOPE -- read this before using the numbers. This script isolates the
# PARAMETER component only: ONE network, ONE allocation, alpha held fixed. It
# says nothing about the allocation/network component, and that omission
# matters, because homogeneous SIR has NO allocation component at all (all
# allocations are exchangeable, so run_fit_array makes a single SIR job) while
# the network has jobs spanning network seeds and allocations. So this cannot
# on its own explain a network AVE that is NARROWER than SIR's -- adding the
# missing component pushes the network the other way.
#
# What it does establish: beta is poorly identified in the network model
# because the outcome barely responds to it (heavy-tailed contact structure
# saturates -- once the well-connected core has burnt through, more beta only
# reaches poorly-connected nodes). AVE is outcome-like and inherits that
# flatness. So a wider beta posterior does NOT imply a wider AVE; the two
# effects very nearly cancel, and comparing raw beta posterior widths across
# models is not informative about AVE width.
#
# Compare two elasticities at each model's own fit, on matched data:
#
#   dlnC1 / dlnbeta   how fast the DATA responds  -> identifiability of beta
#   dlnAVE/ dlnbeta   how fast the AVE responds   -> how much of a wide beta
#                                                    posterior reaches the AVE
#
# Since sd(ln beta) scales roughly as 1/el_C1 at fixed data precision, the
# transfer from data precision to AVE precision goes as el_AVE / el_C1.
#
#   Rscript diag_ave_beta_sensitivity.R

suppressMessages({library(data.table)})
source("utils.R"); source("stoch_model.R"); source("net_sir_events.R")
net_sir_compile()

N <- 800; N_vac <- 400; N_unvac <- N - N_vac
mean_k <- 6; pa <- 1.4
t_star <- 8; init_I <- 2; gamma <- 1; dt <- 0.05
data_C1 <- 160; data_C2 <- 80
I_ini_2g <- c(3, 3)
n_ve_alloc <- 3          # inner allocations, averaged (the estimand we keep)
n_rep_ve   <- 200
cores <- 8

set.seed(1)
cm  <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
adj <- contact_matrix_to_adj(cm)
csr <- adj_to_csr(adj = adj)
set.seed(11); vac <- sample.int(N, N_vac)

# ---- arm outcomes --------------------------------------------------------
arms_net <- function(b, a, nsim = 400, seed = 3) {
  sus <- rep(1, N); sus[vac] <- a
  r <- run_stoch_network(beta = b, N = N, susceptibility = sus, t = t_star,
        vac = vac, gamma = gamma, timepoints = t_star, I_ini = init_I,
        n_sim = nsim, seed = seed, k_mean = mean_k, cores = cores,
        engine = "events", csr = csr)
  setDT(r); f <- r[time == t_star]; c(C1 = mean(f$C1), C2 = mean(f$C2))
}
arms_sir <- function(b, a, nsim = 400, seed = 3) {
  r <- run_stoch_cd_dust(matrix(1, 2, 2), beta = b, N = c(N_unvac, N_vac),
        t = t_star, I_ini = I_ini_2g, susceptibility = c(1, a), gamma = gamma,
        dt = dt, timepoints = t_star, n_sim = nsim, cores = cores, seed = seed)
  setDT(r); f <- r[time == t_star]; c(C1 = mean(f$C1), C2 = mean(f$C2))
}

# ---- AVE, averaging over inner allocations -------------------------------
ave_net <- function(b, a) {
  r <- get_stoch_eate_network(beta = b, susceptibility = c(1, a), f = 0.5, N = N,
        t = t_star, c_ij = cm, adj = adj, k_mean = mean_k, gamma = gamma,
        n_vac = n_ve_alloc, n_rep = n_rep_ve, timepoints = t_star,
        init_I = init_I, mc.cores = cores, inner_cores = 1)
  setDT(r); mean(r[method == "full_stoch", ave])
}
ave_sir <- function(b, a) {
  r <- get_stoch_eate_sir(beta = b, susceptibility = c(1, a), f = 0.5, N = N,
        t = t_star, gamma = gamma, I_ini = I_ini_2g, n_vac = n_ve_alloc,
        n_rep = n_rep_ve, dt = dt, timepoints = t_star, mc.cores = cores)
  setDT(r); mean(r[method == "full_stoch", ave])
}

# ---- calibrate each model to the same data -------------------------------
calibrate <- function(armf, blo, bhi) {
  a <- 0.3; b <- exp((log(blo) + log(bhi)) / 2)
  for (round in 1:3) {
    lo <- log(blo); hi <- log(bhi)
    for (i in 1:15) { m <- (lo + hi) / 2
      if (armf(exp(m), a)[["C1"]] < data_C1) lo <- m else hi <- m }
    b <- exp((lo + hi) / 2)
    alo <- log(0.01); ahi <- log(2)
    for (i in 1:13) { m <- (alo + ahi) / 2
      if (armf(b, exp(m))[["C2"]] < data_C2) alo <- m else ahi <- m }
    a <- exp((alo + ahi) / 2)
  }
  c(beta = b, alpha = a)
}

p_net <- calibrate(arms_net, 0.2, 60)
p_sir <- calibrate(arms_sir, 0.05, 20)
cat(sprintf("network fit: beta = %.3f  alpha = %.4f\n", p_net[["beta"]], p_net[["alpha"]]))
cat(sprintf("SIR     fit: beta = %.3f  alpha = %.4f\n\n", p_sir[["beta"]], p_sir[["alpha"]]))

# ---- sweep beta multiplicatively around each fit --------------------------
mults <- c(0.7, 0.85, 1, 1.18, 1.4)
sweep <- function(lab, p, armf, avef) {
  rbindlist(lapply(mults, function(m) {
    b <- p[["beta"]] * m
    o <- armf(b, p[["alpha"]], 600, 7)
    data.table(model = lab, mult = m, beta = b,
               C1 = o[["C1"]], C2 = o[["C2"]], AVE = avef(b, p[["alpha"]]))
  }))
}
res <- rbind(sweep("network", p_net, arms_net, ave_net),
             sweep("SIR",     p_sir, arms_sir, ave_sir))
fwrite(res, "output/ave_beta_sensitivity.csv")

cat("=== response to beta (alpha held at each model's fit) ===\n")
print(res[, .(model, mult, beta = round(beta, 3), C1 = round(C1, 1),
              C2 = round(C2, 1), AVE = round(AVE, 4))])

# log-log slopes: elasticity, comparable across models and beta scales
el <- res[, {
  fC <- lm(log(C1) ~ log(beta)); fA <- lm(log(AVE) ~ log(beta))
  .(el_C1 = unname(coef(fC)[2]), el_AVE = unname(coef(fA)[2]))
}, by = model]
cat("\n=== elasticities (dln y / dln beta) ===\n")
print(el[, .(model, el_C1 = round(el_C1, 3), el_AVE = round(el_AVE, 3))])

cat("\n  el_C1  small => the data barely constrains beta => WIDE beta posterior\n")
cat("  el_AVE small => a wide beta posterior barely moves AVE => NARROW AVE\n")
cat("\nImplied AVE cv for a beta posterior of relative width sd_lnbeta:\n")
for (s in c(0.1, 0.2, 0.3)) {
  cat(sprintf("  sd(ln beta) = %.2f ->", s))
  for (m in el$model)
    cat(sprintf("  %s: cv(AVE) = %.3f", m, abs(el[model == m, el_AVE]) * s))
  cat("\n")
}
cat("\nWrote output/ave_beta_sensitivity.csv\n")
