# R wrapper for the exact event-driven network SIR (net_sir_events.cpp).
#
# Drop-in alternative to run_stoch_adj() for the NETWORK model: same process,
# simulated exactly (no dt bias) in O(E log V) instead of O(n * max_degree)
# per dt step.
#
# Validity: time-homogeneous rates (scalar beta, constant gamma), waning = 0.
#
#   source("net_sir_events.R")   # compiles on first use

if (!requireNamespace("cpp11", quietly = TRUE))
  stop("net_sir_events.R needs the 'cpp11' package")

# Compile once per session (cheap: a single translation unit).
# cpp_source() injects the generated R wrapper into `env` and returns a
# DLLInfo, so take the function out of the environment rather than off the
# return value.
.net_sir_env <- new.env(parent = globalenv())
net_sir_compile <- function(force = FALSE) {
  if (!force && !is.null(.net_sir_env$net_sir_event_times))
    return(invisible(.net_sir_env$net_sir_event_times))
  cpp11::cpp_source("net_sir_events.cpp", env = .net_sir_env, quiet = TRUE)
  if (is.null(.net_sir_env$net_sir_event_times))
    stop("net_sir_event_times was not registered by cpp11::cpp_source()")
  invisible(.net_sir_env$net_sir_event_times)
}

# Contact matrix (or an existing adjacency list from contact_matrix_to_adj)
# -> 0-based CSR in the TRANSMIT direction.
#
# Direction matters: get_conact_matrix_pl returns an ASYMMETRIC (directed)
# matrix. In stoch_mod_adj.R the FOI on node i sums over neighbors[, i], i.e.
# row i lists who infects i -- a RECEIVE list. The event-driven simulator
# walks outward from an infectious node, so it needs the transpose: who node k
# can infect, = {i : k in recv(i)} = column k. Using the receive list to
# transmit silently simulates the reversed graph.
adj_to_csr <- function(contact_matrix = NULL, adj = NULL) {
  if (is.null(adj)) {
    n    <- nrow(contact_matrix)
    recv <- lapply(seq_len(n), function(i) which(contact_matrix[i, ] != 0))
  } else {
    n    <- ncol(adj$neighbors)
    recv <- lapply(seq_len(n), function(i) {
      k <- which(adj$mask[, i] == 1L)
      if (!length(k)) integer(0) else as.integer(adj$neighbors[k, i])
    })
  }
  # Transpose recv -> transmit: edge (k infects i) for each i, k in recv[[i]].
  to    <- rep.int(seq_len(n), lengths(recv))          # i (receiver)
  from  <- as.integer(unlist(recv, use.names = FALSE)) # k (transmitter)
  ord   <- order(from, method = "radix")
  deg   <- tabulate(from, nbins = n)
  list(n   = n,
       nbr = as.integer(to[ord]) - 1L,                 # 0-based receivers
       ptr = as.integer(c(0L, cumsum(deg))))
}

# Returns a data.table matching run_stoch_network()'s output: time, sim,
# C1 (cumulative unvaccinated cases), C2 (cumulative vaccinated cases).
#
# `beta` is on the SAME scale run_stoch_adj receives, i.e. run_stoch_network
# has already multiplied the user's beta by N / k_mean.
run_stoch_network_events <- function(beta, N, susceptibility = c(1, 1),
                                     t, vac, csr, gamma = 1,
                                     timepoints = NULL, I_ini = 2,
                                     transmissibility = NULL,
                                     n_sim = 100, seed = 1L, k_mean = 6) {
  net_sir_compile()
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n <- csr$n
  sus <- rep(susceptibility[1], n); sus[vac] <- susceptibility[2]
  tr  <- if (is.null(transmissibility)) rep(1, n) else transmissibility
  seeds <- seq_len(min(I_ini, n)) - 1L        # matches run_stoch_network's I_ini

  inf <- .net_sir_env$net_sir_event_times(
    n = as.integer(n), nbr = csr$nbr, ptr = csr$ptr,
    susceptibility = as.numeric(sus), transmissibility = as.numeric(tr),
    beta = as.numeric(N * beta / k_mean), gamma = as.numeric(gamma),
    t_max = as.numeric(max(timepoints)),
    seeds = as.integer(seeds), n_sim = as.integer(n_sim),
    seed = as.integer(seed))

  is_vac <- logical(n); is_vac[vac] <- TRUE
  data.table::rbindlist(lapply(timepoints, function(tp) {
    hit <- inf <= tp                          # [n_sim, n] cumulative incidence
    data.table::data.table(
      time = tp, sim = seq_len(n_sim),
      C1 = as.numeric(rowSums(hit[, !is_vac, drop = FALSE])),
      C2 = as.numeric(rowSums(hit[,  is_vac, drop = FALSE])))
  }))
}
