library(dplyr)
library(data.table)

# ---------------------------------------------------------------------------
# Frailty distribution
# ---------------------------------------------------------------------------

get_frailty <- function(mean=0.5, sd=1, n=100){
  variance <- sd^2
  common_factor <- (mean * (1 - mean) / variance) - 1
  shape1 <- mean * common_factor
  shape2 <- (1 - mean) * common_factor
  x <- (seq(0, 1, length.out=n+1) + 1/(n)/2)[1:n]
  p <- dbeta(x, shape1, shape2)
  return(list(x=x, p=p/sum(p)))
}

# ---------------------------------------------------------------------------
# Transmission rate from R0
# ---------------------------------------------------------------------------

get_beta <- function(R, alpha, sd, sd_trans=0, f=0.5, N=1000, n_frailty=100, gamma=1/2){
  # NGM[i,j] = mm[i,j] * sus[i] * S_ini[i] * tr[j] / (gamma * N_total).
  # With sd_trans > 0 the trans frailty is rank-correlated with the sus
  # frailty exactly as in run_stoch_frailty_cd, so the calibration matches
  # the simulator.
  sd_pop <- if (sd > 0) sd else sd_trans
  if (sd_pop > 0) {
    fr        <- get_frailty(sd=sd_pop, n=n_frailty)
    # Normalised to population-mean 1 so beta=R0/gamma in the homogeneous
    # limit (matches run_stoch_frailty_cd / get_frailty_eate).
    sus_bin   <- if (sd > 0) {
      raw <- exp(2.5 * fr$x); raw / sum(fr$p * raw)
    } else rep(1, n_frailty)
    trans_bin <- if (sd_trans > 0) {
      raw <- if (sd_trans == sd_pop) {
        exp(2.5 * fr$x)
      } else {
        cf_pop <- (0.25 / sd_pop^2)   - 1
        cf_t   <- (0.25 / sd_trans^2) - 1
        ranks  <- pbeta(fr$x, 0.5 * cf_pop, 0.5 * cf_pop)
        exp(2.5 * qbeta(ranks, 0.5 * cf_t, 0.5 * cf_t))
      }
      raw / sum(fr$p * raw)
    } else {
      rep(1, n_frailty)
    }
    S_ini <- c(2*N*(1-f) * fr$p, 2*N*f * fr$p)
    sus   <- c(sus_bin,   alpha * sus_bin)
    tr    <- c(trans_bin, trans_bin)              # vaccine acts on sus only
  } else {
    S_ini <- c(2*N*(1-f), 2*N*f)
    sus   <- c(1, alpha)
    tr    <- c(1, 1)
  }
  n_g <- length(S_ini)
  # Unit mixing entries (NOT /n_g) to match run_stoch_frailty_cd /
  # run_frailty_cd convention; combined with frailty arrays normalised to
  # mean 1, this gives R0 = beta/gamma in the homogeneous limit.
  mm  <- matrix(1, nrow=n_g, ncol=n_g)
  ng  <- outer(sus * S_ini, tr) * mm / (gamma * 2 * N)
  eig <- Re(eigen(ng, only.values=TRUE)$values[1])
  R / eig
}

# ---------------------------------------------------------------------------
# Next-generation matrix
# ---------------------------------------------------------------------------

cij_NGM <- function(c_ij, N, susceptibility, transmisibility, gamma=1/3, norm_contacts=NULL){
  if(is.null(norm_contacts)){
    norm <- c_ij %*% N/sum(N)
  }else{
    N_conts <- as.numeric(norm_contacts %*% N)
    norm <- c_ij %*% N/sum(N)/N_conts*sum(N_conts)
  }
  c_ij <- c_ij/as.numeric(norm)
  NGM <- c_ij %*% diag(transmisibility)*N/sum(N)*susceptibility
  beta_R <- Re(eigen(NGM, only.values=T)$values[1]/gamma)
  return(list(c_ij=c_ij, NGM=NGM, beta_R=beta_R))
}

# ---------------------------------------------------------------------------
# Pareto (power-law) contact network
# ---------------------------------------------------------------------------

get_conact_matrix_pl <- function(N, alpha, mean_k=6){
  propensities <- Pareto::rPareto(N, alpha=alpha, t=1)
  propensities <- propensities / mean(propensities)
  contacts <- matrix(0, nrow=N, ncol=N)
  for(i in 1:N){
    for(j in 1:N){
      contacts[i,j] <- rbinom(1, 1, min(mean_k*propensities[i]*propensities[j]/N, 1))
    }
  }
  return(contacts)
}

# ---------------------------------------------------------------------------
# Adjacency list from binary contact matrix (for sparse odin2 model)
# ---------------------------------------------------------------------------

contact_matrix_to_adj <- function(contact_matrix) {
  n <- nrow(contact_matrix)
  adj <- lapply(1:n, function(i) which(contact_matrix[i, ] != 0))
  degree <- lengths(adj)
  max_degree <- max(degree, 1L)
  neighbors <- matrix(1L,  nrow = max_degree, ncol = n)
  mask      <- matrix(0L,  nrow = max_degree, ncol = n)
  for (i in seq_len(n)) {
    if (degree[i] > 0) {
      neighbors[1:degree[i], i] <- adj[[i]]
      mask[1:degree[i], i]      <- 1L
    }
  }
  list(neighbors = neighbors, mask = mask, max_degree = max_degree)
}

# Weighted variant: returns the actual contact_matrix[i, j] entries in a
# weights matrix [max_degree, n] alongside neighbors. Used by stoch_ind.R.
contact_matrix_to_weighted_adj <- function(contact_matrix) {
  n <- nrow(contact_matrix)
  adj <- lapply(1:n, function(i) which(contact_matrix[i, ] != 0))
  degree <- lengths(adj)
  max_degree <- max(degree, 1L)
  neighbors <- matrix(1L, nrow = max_degree, ncol = n)
  weights   <- matrix(0,  nrow = max_degree, ncol = n)
  for (i in seq_len(n)) {
    if (degree[i] > 0) {
      neighbors[1:degree[i], i] <- adj[[i]]
      weights[1:degree[i], i]   <- contact_matrix[i, adj[[i]]]
    }
  }
  list(neighbors = neighbors, weights = weights, max_degree = max_degree)
}

# ---------------------------------------------------------------------------
# Plot theme
# ---------------------------------------------------------------------------

add_theme <- function(q){
  q + theme_bw() + theme(text = element_text(size=8)) + scale_size_identity()
}


# ---------------------------------------------------------------------------
# Approximate posterior covariance (J^-1 Sigma J^-T)
# ---------------------------------------------------------------------------
#
# Given a stochastic simulator and a fitted (beta, alpha), estimate the
# Gaussian/Laplace posterior covariance via the sandwich
#     Cov(theta_hat) ~ J^-1 Sigma J^-T
# where
#   J[i, j] = d E[C_i] / d theta_j  (central finite differences),
#   Sigma   = Cov(C_i) / n_sim      (the standard error of mean estimates),
# computed with common random numbers (CRN) — the same dust2 seed is
# used for all 5 simulator calls so the finite-difference noise drops
# from O(1/sqrt(n_sim)) to O(h^2). The base run also supplies Sigma so
# we don't have to pool the ± perturbed runs (which would inflate Sigma
# with a between-group component).
#
# `simulator(beta, alpha, n_sim, seed)` must return a data.frame with
# columns C1, C2 (one row per replicate at the focal time).
#
# Returns a list with
#   cov   2x2 posterior covariance of (beta, alpha)
#   J     2x2 Jacobian (rows = (C1, C2), cols = (beta, alpha))
#   Sigma 2x2 covariance of the means of (C1, C2)
#   base  the centre-point simulator output (for diagnostics)
#   sd    named vector (beta, alpha) of marginal posterior SDs
estimate_posterior_cov <- function(simulator, beta, alpha,
                                   n_sim = 1000, seed = 42L,
                                   h_rel = 0.01) {
  hb <- beta  * h_rel
  ha <- alpha * h_rel

  # Boundary safety: if symmetric central step would leave the support,
  # shrink h. (Both beta and alpha are required positive by the simulators.)
  if (beta - hb <= 0)  hb <- beta  / 2
  if (alpha - ha <= 0) ha <- alpha / 2

  base     <- simulator(beta,           alpha,            n_sim, seed = seed)
  pb_plus  <- simulator(beta + hb,      alpha,            n_sim, seed = seed)
  pb_minus <- simulator(beta - hb,      alpha,            n_sim, seed = seed)
  pa_plus  <- simulator(beta,           alpha + ha,       n_sim, seed = seed)
  pa_minus <- simulator(beta,           alpha - ha,       n_sim, seed = seed)

  J <- matrix(
    c((mean(pb_plus$C1) - mean(pb_minus$C1)) / (2 * hb),
      (mean(pb_plus$C2) - mean(pb_minus$C2)) / (2 * hb),
      (mean(pa_plus$C1) - mean(pa_minus$C1)) / (2 * ha),
      (mean(pa_plus$C2) - mean(pa_minus$C2)) / (2 * ha)),
    nrow = 2, byrow = FALSE,
    dimnames = list(c("C1", "C2"), c("beta", "alpha"))
  )

  Sigma <- cov(cbind(base$C1, base$C2)) / n_sim
  dimnames(Sigma) <- list(c("C1", "C2"), c("C1", "C2"))

  Jinv <- tryCatch(solve(J),
                   error = function(e) {
                     warning("Jacobian is singular; returning NA covariance. ",
                             "Inspect $J for the cause.")
                     matrix(NA_real_, 2, 2)
                   })
  cov_par <- Jinv %*% Sigma %*% t(Jinv)
  dimnames(cov_par) <- list(c("beta", "alpha"), c("beta", "alpha"))

  list(cov   = cov_par,
       J     = J,
       Sigma = Sigma,
       base  = base,
       sd    = sqrt(c(beta = cov_par[1, 1], alpha = cov_par[2, 2])))
}
