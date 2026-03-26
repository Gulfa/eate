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

get_beta <- function(R, alpha, sd, f=0.5, N=1000, n_frailty=100, gamma=1/2){
  fr <- get_frailty(sd=sd, n=n_frailty)
  frailty <- exp(2.5*fr$x)
  params <- list(
    n=2*n_frailty,
    S_ini=c(2*N*(1-f)*fr$p, 2*N*f*fr$p),
    susceptibility=c(frailty, alpha*frailty))
  mm <- matrix(1, nrow=params$n, ncol=params$n)/params$n
  ng <- sweep(mm, MARGIN=1, params$susceptibility*params$S_ini/gamma/(2*N), `*`)
  eig <- Re(eigen(ng, only.values=T)$values[1])
  beta <- R / eig
  return(beta)
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

# ---------------------------------------------------------------------------
# Plot theme
# ---------------------------------------------------------------------------

add_theme <- function(q){
  q + theme_bw() + theme(text = element_text(size=8)) + scale_size_identity()
}
