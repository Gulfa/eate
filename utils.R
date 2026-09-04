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


# Sparse directed Pareto (Chung-Lu style) network sampler. Draws exactly
# the same Bernoulli edges as get_conact_matrix_pl but never materialises
# the dense N x N matrix.
#
# Uses Miller-Hagberg-style rejection sampling with a *row-fixed*
# envelope: for source i, sort target propensities descending and use
#     p_max = min(c * s[i] * s_sort[1], 1)
# as the acceptance rate of the outer Bernoulli process, then accept
# each candidate position pos with probability p_true(pos) / p_max
# where p_true(pos) = c * s[i] * s_sort[pos]. Because s_sort[pos] is
# monotonically decreasing, p_true <= p_max is guaranteed and every
# candidate integrates to the right marginal Bernoulli probability.
#
# Vectorised across positions: for row i we pre-sample enough Geom(p_max)
# gaps at once (with a safety multiplier), take a cumsum to get all
# candidate positions in one shot, then vectorised runif() gives all
# accept/reject decisions. Row work is O(p_max * N + #accepts), and
# summed over rows the total work is O(mean_k * max(s) * N) — well
# below O(N^2) for moderate pl_alpha.
#
# Returns the same list shape as contact_matrix_to_adj:
#   neighbors [max_degree, N] integer, mask [max_degree, N] integer,
#   max_degree scalar, plus propensities (length N) and degree (length N).
#
# NOTE: an earlier attempt at MH with mid-row envelope refresh had a
# subtle bias in the directed accept ratio and gave systematically wrong
# CIR curves — see the fix commit. This version uses the row-fixed p_max
# variant from the MH 2011 paper, which is provably unbiased.
# symmetric = TRUE (default) draws each unordered pair ONCE and mirrors it, so
# the contact graph is undirected: if i can infect j then j can infect i. The
# original code drew every row independently, so p_ij and p_ji were separate
# Bernoulli trials of the same probability and agreed only by chance --
# reciprocity came out at ~1.3%, i.e. an almost fully one-directional graph.
# That silently suppresses any feedback effect (measured: the indirect effect
# of vaccinating i back onto i was +0.0001 directed vs -0.0254 symmetric).
#
# The construction is the standard Chung-Lu one: sample only the upper triangle
# in sorted order and mirror. Same skip-sampling, same O(N + m) cost, same
# expected mean_k and power-law degrees -- unlike post-hoc symmetrisation with
# max(A, A') which doubles the degree, or A & A' which collapses it.
sample_pareto_adj <- function(N, alpha, mean_k = 6, seed = NULL,
                              allow_self_loops = TRUE,
                              oversample = 1.5, symmetric = TRUE) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
  }
  s   <- Pareto::rPareto(N, alpha = alpha, t = 1)
  s   <- s / mean(s)
  c_const <- mean_k / N

  # Sort target propensities descending so envelope is non-increasing.
  ord    <- order(s, decreasing = TRUE)
  s_sort <- s[ord]
  max_s  <- s_sort[1]

  neigh <- vector("list", N)
  # Symmetric: walk in SORTED order and only consider partners of higher rank,
  # so each unordered pair is offered exactly once. `start` is that rank floor
  # (0 in the directed case, where every row scans all N candidates).
  rank_seq <- if (symmetric) seq_len(N) else seq_len(N)
  for (r in rank_seq) {
    i     <- if (symmetric) ord[r] else r
    start <- if (symmetric) r else 0L
    si    <- s[i]
    p_max <- min(c_const * si * max_s, 1)

    if (p_max <= 0) { neigh[[i]] <- integer(0L); next }

    if (p_max >= 1) {
      # Envelope saturated — no advantage to skipping. Fall back to
      # direct Bernoulli. Only affects the top few rows for heavy tails.
      if (symmetric) {
        cand  <- if (start < N) (start + 1L):N else integer(0L)
        p_vec <- pmin(c_const * si * s_sort[cand], 1)
        out   <- if (length(cand)) ord[cand[rbinom(length(cand), 1L, p_vec) == 1L]]
                 else integer(0L)
      } else {
        p_vec <- pmin(c_const * si * s, 1)
        out   <- which(rbinom(N, 1L, p_vec) == 1L)
      }
    } else {
      # Vectorised Geom(p_max) skip + cumsum to enumerate candidate
      # positions in sorted order. Pre-allocate `oversample * expected`
      # gaps; refill loop kicks in the rare cases we underestimate.
      pos_out <- integer(0L)
      last_pos <- start
      repeat {
        n_batch <- max(64L, ceiling(oversample * p_max * (N - last_pos)))
        gaps    <- rgeom(n_batch, p_max)
        pos_vec <- last_pos + cumsum(gaps + 1L)
        keep    <- pos_vec <= N
        pos_out <- c(pos_out, pos_vec[keep])
        # If any candidate exceeded N, we've walked past — done.
        if (!all(keep)) break
        last_pos <- pos_vec[length(pos_vec)]
      }
      # Rejection step at the true (smaller) probability.
      if (length(pos_out) > 0L) {
        p_true <- c_const * si * s_sort[pos_out]
        accept <- runif(length(pos_out)) < p_true / p_max
        out    <- ord[pos_out[accept]]
      } else {
        out <- integer(0L)
      }
    }

    if (!allow_self_loops) out <- out[out != i]
    if (symmetric) {
      # Mirror: record the edge on both endpoints.
      out <- out[out != i]
      neigh[[i]] <- c(neigh[[i]], out)
      for (j in out) neigh[[j]] <- c(neigh[[j]], i)
    } else {
      neigh[[i]] <- out
    }
  }

  if (!is.null(seed)) {
    if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv)
    else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }

  degree     <- lengths(neigh)
  max_degree <- max(degree, 1L)
  neighbors  <- matrix(1L, nrow = max_degree, ncol = N)
  mask       <- matrix(0L, nrow = max_degree, ncol = N)
  for (i in seq_len(N)) {
    if (degree[i] > 0L) {
      neighbors[seq_len(degree[i]), i] <- neigh[[i]]
      mask[     seq_len(degree[i]), i] <- 1L
    }
  }
  list(neighbors    = neighbors,
       mask         = mask,
       max_degree   = as.integer(max_degree),
       propensities = s,
       degree       = degree)
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

  # Variance of a SINGLE realisation of (C1, C2) — the residual variance
  # that appears in the sandwich estimator for a moment-matching MLE.
  # The n_sim replicates are used to *estimate* this variance from an
  # ensemble; each row of `base` is one realisation, so cov() already
  # gives Sigma_1. Do NOT divide by n_sim — that would give Var(mean of
  # n_sim draws) instead, understating posterior SDs by sqrt(n_sim) ~ 30x
  # for the defaults (n_sim = 1000). The fit is against ONE observed
  # (data_C1, data_C2), so Sigma_1 is what the sandwich wants.
  Sigma <- cov(cbind(base$C1, base$C2))
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
