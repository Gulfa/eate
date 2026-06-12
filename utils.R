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
    sus_bin   <- if (sd > 0) exp(2.5 * fr$x) else rep(1, n_frailty)
    trans_bin <- if (sd_trans > 0) {
      if (sd_trans == sd_pop) {
        exp(2.5 * fr$x)
      } else {
        cf_pop <- (0.25 / sd_pop^2)   - 1
        cf_t   <- (0.25 / sd_trans^2) - 1
        ranks  <- pbeta(fr$x, 0.5 * cf_pop, 0.5 * cf_pop)
        exp(2.5 * qbeta(ranks, 0.5 * cf_t, 0.5 * cf_t))
      }
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
  mm  <- matrix(1, nrow=n_g, ncol=n_g) / n_g
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
