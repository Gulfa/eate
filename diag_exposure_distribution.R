# Why the network needs a much smaller alpha than SIR for the same trial data
# ---------------------------------------------------------------------------
# Fitted to identical data (AR 0.40 unvaccinated / 0.20 vaccinated at t* = 8),
# the homogeneous SIR wants alpha = 0.415 and the pl_alpha = 1.4 network wants
# alpha = 0.240 -- a 1.7x stronger vaccine asserted from the same trial.
#
# The reason is the shape of the exposure distribution, and it is invisible in
# the marginal per-node infection probability. On a network a LOW attack rate
# does not mean everyone has low exposure: it means most people are never
# reached at all while the few the epidemic does reach are heavily exposed.
# Averaging over realisations first smears those two populations into a
# moderate-looking middle value that describes nobody.
#
# This extracts the per-(realisation, node) cumulative exposure Lambda and
# solves for the alpha that produces CIR = 0.5 under three descriptions of it:
# the correct per-realisation one, the misleading marginal one, and a
# homogeneous population at the same attack rate.
#
# Consequence for the paper: the fitted alpha is not a transportable
# "biological" parameter. It absorbs the exposure distribution of whatever
# model produced it. It also means 1 - alpha is NOT an attainable ceiling on a
# network -- lowering the attack rate removes people from exposure rather than
# reducing everyone's, so the exposed stay at Lambda ~ O(1) and the ratio never
# approaches alpha.
#
#   Rscript diag_exposure_distribution.R
# ---------------------------------------------------------------------------

suppressMessages({library(data.table)})
source("utils.R"); source("stoch_model.R")

N <- 800; mean_k <- 6; pa <- 1.4; t_star <- 8; gamma <- 1; dt <- 0.1
beta <- 3.536                       # the network's own calibration to the data
tps <- seq(1, t_star, 1); n_t <- length(tps); n_rep <- 250

set.seed(1)
c_ij <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
adj  <- contact_matrix_to_adj(c_ij)

# Everyone unvaccinated: the bare exposure field the vaccine would act against.
raw <- run_stoch_adj(c_ij, beta = N * beta / mean_k, t = t_star,
                     I_ini = c(rep(1L, 2), rep(0L, N - 2)),
                     susceptibility = rep(1, N), gamma = gamma, dt = dt,
                     timepoints = tps, n_sim = n_rep, cores = 8, adj = adj)
setDT(raw)

I_mat <- array(0, dim = c(n_t, n_rep, N))
for (k in seq_len(N))
  I_mat[, , k] <- .dt_col_to_t_rep_matrix(raw[[paste0("I", k)]], n_t, n_rep)
cf <- array(0, dim = c(n_t, n_rep, N))
for (r in seq_len(n_rep))
  cf[, r, ] <- .cum_trapz((I_mat[, r, ] %*% t(c_ij)) * (beta / mean_k), tps)

L  <- as.vector(cf[n_t, , ])                       # per (realisation, node)
Pi <- colMeans(1 - exp(-cf[n_t, , ]))              # marginal per node
Lm <- -log(pmax(1 - Pi, 1e-12))
AR <- mean(1 - exp(-L))
Lhom <- -log(1 - AR)

cat(sprintf("network pl_alpha = %.1f, N = %d, beta = %.3f, t* = %d\n\n",
            pa, N, beta, t_star))
cat(sprintf("per-realisation exposure (%d node x realisation pairs)\n", length(L)))
cat(sprintf("  never reached (L < 0.01): %.1f%%      L > 1: %.1f%%\n",
            100 * mean(L < 0.01), 100 * mean(L > 1)))
cat(sprintf("  among the exposed: mean %.2f  median %.2f  q90 %.2f\n",
            mean(L[L > 0.01]), median(L[L > 0.01]), quantile(L[L > 0.01], .9)))
cat(sprintf("  overall attack rate %.3f\n", AR))
cat(sprintf("\nmarginal Lambda per node: mean %.3f  median %.3f  max %.2f",
            mean(Lm), median(Lm), max(Lm)))
cat("   <- tame, and wrong\n")

solve_a <- function(x, target = 0.5)
  uniroot(function(a) sum(1 - exp(-a * x)) / sum(1 - exp(-x)) - target,
          c(1e-6, 1))$root

cat("\nalpha required to produce CIR = 0.5:\n")
cat(sprintf("  per-realisation exposure (correct)    %.4f   vs network fit 0.2403\n",
            solve_a(L)))
cat(sprintf("  homogeneous at the same attack rate   %.4f   vs SIR     fit 0.4147\n",
            solve_a(rep(Lhom, length(L)))))
cat(sprintf("  marginal Lambda per node (misleading) %.4f   explains neither\n",
            solve_a(Lm)))
