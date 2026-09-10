# Is the frozen-field approximation self-consistent on a network?
# ---------------------------------------------------------------------------
# get_stoch_eate_network forms a HYBRID contrast: factually-matching individuals
# contribute their empirical infection probability P_fac, flipped individuals
# contribute a frozen counterfactual built from the per-node cumulative force of
# infection Lambda_i,
#     P_unvac_cf = 1 - E_r[exp(-Lambda_i)]      P_vac_cf = 1 - E_r[exp(-alpha*Lambda_i)]
#     num   = sum_vac P_fac + sum_unvac P_vac_cf
#     denom = sum_unvac P_fac + sum_vac P_unvac_cf
#
# That is only coherent if the frozen field reproduces the factual world for the
# people who were NOT flipped, i.e. if
#     P_fac[i] ~ P_unvac_cf[i]  for unvaccinated i
#     P_fac[i] ~ P_vac_cf[i]    for vaccinated i
# On a network the exposures behind Lambda_i are strongly correlated (a few
# infected neighbours, not an independent aggregate), so 1 - exp(-Lambda) may
# not be the right first-passage probability at all. If it is not, the hybrid
# mixes two inconsistent quantities and its ratio is biased.
#
# This rebuilds the internals of get_stoch_eate_network and reports, per arm,
# the empirical vs frozen-field probabilities, plus the hybrid EATE against a
# PURE counterfactual EATE (both arms from the frozen field) which is internally
# consistent by construction. Run at two coverages, since the question that
# raised this was why network VE barely moves with coverage.
#
#   Rscript diag_frozen_field.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table)})
source("utils.R"); source("stoch_model.R")

N <- 800; pa <- 1.4; mean_k <- 6
t_star <- 8; init_I <- 2; gamma <- 1; dt <- 0.1
beta <- 4.90; alpha <- 0.19          # ~the pa=1.4 fits from diag_ave_refit
n_rep <- 300
timepoints <- seq(1, t_star, 1)
n_t <- length(timepoints)
covs <- c(0.25, 0.75)

set.seed(1)
c_ij <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
adj  <- contact_matrix_to_adj(c_ij)

one <- function(f, seed = 11) {
  set.seed(seed)
  vac <- sample(seq_len(N), round(f * N))
  non_vac <- setdiff(seq_len(N), vac)

  susept <- rep(1, N); susept[vac] <- alpha
  I_ini_vec <- c(rep(1L, init_I), rep(0L, N - init_I))
  raw <- run_stoch_adj(c_ij, beta = N * beta / mean_k, t = t_star,
                       I_ini = I_ini_vec, susceptibility = susept,
                       gamma = gamma, dt = dt, timepoints = timepoints,
                       n_sim = n_rep, cores = 8, adj = adj)
  setDT(raw)

  P_factual <- matrix(0, n_t, N)
  I_mat     <- array(0, dim = c(n_t, n_rep, N))
  for (k in seq_len(N)) {
    P_factual[, k] <- rowMeans(.dt_col_to_t_rep_matrix(raw[[paste0("C", k)]],
                                                       n_t, n_rep))
    I_mat[, , k]   <- .dt_col_to_t_rep_matrix(raw[[paste0("I", k)]], n_t, n_rep)
  }
  cum_foi <- array(0, dim = c(n_t, n_rep, N))
  for (r in seq_len(n_rep)) {
    FI_r <- (I_mat[, r, ] %*% t(c_ij)) * (beta / mean_k)
    cum_foi[, r, ] <- .cum_trapz(FI_r, timepoints)
  }

  it   <- n_t                                   # read at t_star
  cfi  <- cum_foi[it, , ]                       # [n_rep, N]
  P_vac_cf   <- 1 - colMeans(exp(-alpha * cfi))
  P_unvac_cf <- 1 - colMeans(exp(-cfi))
  P_fac      <- P_factual[it, ]

  # hybrid (what the code does) vs pure counterfactual (internally consistent)
  num_h <- sum(P_fac[vac])     + sum(P_vac_cf[non_vac])
  den_h <- sum(P_fac[non_vac]) + sum(P_unvac_cf[vac])
  num_c <- sum(P_vac_cf);  den_c <- sum(P_unvac_cf)

  list(f = f,
       fac_unvac = mean(P_fac[non_vac]),   ff_unvac = mean(P_unvac_cf[non_vac]),
       fac_vac   = mean(P_fac[vac]),       ff_vac   = mean(P_vac_cf[vac]),
       lam_mean  = mean(colMeans(cfi)),
       VE_hybrid = 1 - num_h / den_h,
       VE_pure   = 1 - num_c / den_c,
       VE_fac    = 1 - (mean(P_fac[vac]) / mean(P_fac[non_vac])))
}

res <- rbindlist(lapply(covs, function(f) as.data.table(one(f))))

cat(sprintf("network pl_alpha=%.1f  N=%d  beta=%.2f  alpha=%.2f  t*=%d  n_rep=%d\n\n",
            pa, N, beta, alpha, t_star, n_rep))
cat("=== is the frozen field self-consistent for the NON-flipped people? ===\n")
cat(sprintf("%-6s %-22s %-22s %s\n", "cov",
            "unvaccinated arm", "vaccinated arm", "mean Lambda"))
cat(sprintf("%-6s %-10s %-10s %-10s %-10s\n", "", "empirical", "frozen",
            "empirical", "frozen"))
for (i in seq_len(nrow(res))) with(res[i], {
  cat(sprintf("%-6.2f %-10.4f %-10.4f %-10.4f %-10.4f %.3f\n",
              f, fac_unvac, ff_unvac, fac_vac, ff_vac, lam_mean))
})

cat("\n=== what that does to the estimand ===\n")
print(res[, .(coverage = f,
              VE_hybrid = round(VE_hybrid, 4),
              VE_pure   = round(VE_pure, 4),
              VE_factual_ratio = round(VE_fac, 4))])
cat(sprintf("\nspan across coverage:  hybrid %.4f   pure counterfactual %.4f\n",
            diff(range(res$VE_hybrid)), diff(range(res$VE_pure))))
cat(sprintf("1 - alpha = %.3f\n", 1 - alpha))
cat("\nIf empirical and frozen disagree, the hybrid is mixing inconsistent\n",
    "quantities and VE_hybrid should not be read as the causal contrast.\n")
fwrite(res, "output/frozen_field_check.csv")
cat("\nWrote output/frozen_field_check.csv\n")
