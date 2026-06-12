library(adaptivetau)
library(data.table)

stoch_model_cd     <- odin2::odin("stoch_sir.R")
stoch_model_adj    <- odin2::odin("stoch_mod_adj.R")
stoch_model_linear <- odin2::odin("stoch_linear.R")
stoch_model_ind    <- odin2::odin("stoch_ind.R")

transitions_sir <- list(
  c(S1=-1, I1=1, C1=1), c(I1=-1, R1=1),
  c(S2=-1, I2=1, C2=1), c(I2=-1, R2=1)
)

transitions_linear <- list(
  c(S1=-1, C1=1),
  c(S2=-1, C2=1)
)

sir_rates <- function(x, p, t) {
  I_total <- x["I1"] + x["I2"]
  c(p$beta * x["S1"] * I_total / p$N,
    p$gamma * x["I1"],
    p$alpha * p$beta * x["S2"] * I_total / p$N,
    p$gamma * x["I2"])
}

linear_rates <- function(x, p, t) {
  c(p$beta * x["S1"],
    p$alpha * p$beta * x["S2"])
}

run_sir <- function(N_cont, N_vac, beta, alpha, gamma=1/7, I0_cont=1, I0_vac=0, t=100,
                    timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  params <- list(N=N_cont + N_vac, beta=beta, alpha=alpha, gamma=gamma)
  init <- c(S1=N_cont - I0_cont, I1=I0_cont, R1=0, C1=0,
            S2=N_vac  - I0_vac,  I2=I0_vac,  R2=0, C2=0)

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions_sir, sir_rates, params, tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

run_linear <- function(N_cont, N_vac, beta, alpha, t=100,
                       timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  params <- list(beta=beta, alpha=alpha)
  init <- c(S1=N_cont, C1=0, S2=N_vac, C2=0)

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions_linear, linear_rates, params, tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

run_stoch_cd_ctmc <- function(mixing_matrix, beta, N, t, I_ini,
                              susceptibility=NULL, transmissibility=NULL,
                              gamma=1/3, timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n <- nrow(mixing_matrix)
  if (is.null(susceptibility))   susceptibility   <- rep(1, n)
  if (is.null(transmissibility)) transmissibility <- rep(1, n)
  N_total <- sum(N)

  S_names <- paste0("S", 1:n)
  I_names <- paste0("I", 1:n)
  R_names <- paste0("R", 1:n)
  C_names <- paste0("C", 1:n)

  transitions <- unlist(recursive=FALSE, lapply(1:n, function(i) {
    inf <- setNames(c(-1L,  1L, 1L), c(S_names[i], I_names[i], C_names[i]))
    rec <- setNames(c(-1L,  1L),     c(I_names[i], R_names[i]))
    list(inf, rec)
  }))

  rate_func <- local({
    mm  <- mixing_matrix
    sus <- susceptibility
    tr  <- transmissibility
    function(x, p, t) {
      I_vec <- x[I_names]
      foi   <- as.numeric(mm %*% (I_vec * tr)) / N_total
      rates <- numeric(2 * n)
      rates[seq(1, 2*n, 2)] <- pmax(0, sus * beta * x[S_names] * foi)
      rates[seq(2, 2*n, 2)] <- pmax(0, gamma * I_vec)
      rates
    }
  })

  init <- setNames(
    c(rbind(N - I_ini, I_ini, rep(0, n), rep(0, n))),
    c(rbind(S_names,   I_names, R_names,   C_names))
  )

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions, rate_func, list(), tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

# Discrete-time tau-leaping equivalent of run_stoch_cd_ctmc, backed by the
# odin2/dust2 model in stoch_sir.R. The R-level signature mirrors
# run_stoch_cd_ctmc; the odin model's parameter interface mirrors det_mod_cd.R
# (beta_day, waning, ...). Replication uses dust2's native n_particles +
# n_threads (one shared system, no R-level forking).
#
# `dt` is the tau-leaping step size. The default 0.1 gives results within a
# few percent of run_stoch_cd_ctmc for typical parameters. With dt close to
# 1/gamma (e.g. dt=1 and gamma~1) expect noticeable bias on cumulative
# infections — shrink dt for higher rates.
run_stoch_cd_dust <- function(mixing_matrix, beta, N, t, I_ini,
                              susceptibility=NULL, transmissibility=NULL,
                              gamma=1/3, waning=0, dt=0.1,
                              timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n <- nrow(mixing_matrix)
  if (is.null(susceptibility))   susceptibility   <- rep(1, n)
  if (is.null(transmissibility)) transmissibility <- rep(1, n)

  # One beta value per dt step. Total steps to reach the last requested time.
  n_steps_total <- as.integer(ceiling(max(timepoints) / dt))
  beta_day_vec  <- if (length(beta) == 1L) rep(beta, n_steps_total) else beta
  N_steps       <- length(beta_day_vec)

  params <- list(
    n               = as.integer(n),
    S_ini           = N - I_ini,
    I_ini           = I_ini,
    mixing_matrix   = mixing_matrix,
    beta_day        = beta_day_vec,
    susceptibility  = susceptibility,
    transmisibility = transmissibility,
    N_steps         = as.integer(N_steps),
    waning          = waning,
    gamma           = gamma
  )

  sys <- dust2::dust_system_create(stoch_model_cd, params,
                                   n_particles = n_sim,
                                   n_threads   = cores,
                                   dt          = dt,
                                   time        = 0)
  dust2::dust_system_set_state_initial(sys)
  raw <- dust2::dust_system_simulate(sys, timepoints)
  # raw shape: [n_states, n_particles, n_times] (or [n_states, n_times] if n_sim==1)
  if (length(dim(raw)) == 2L) {
    raw <- array(raw, dim = c(dim(raw)[1L], 1L, dim(raw)[2L]))
  }
  # State layout (declaration order in stoch_sir.R):
  #   S[1..n], I[1..n], R[1..n], C[1..n], exposure[1..n]
  n_t <- length(timepoints)
  out <- data.table(
    time = rep(timepoints,     each  = n_sim),
    sim  = rep(seq_len(n_sim), times = n_t)
  )
  offsets <- c(S=0L, I=n, R=2L*n, C=3L*n)
  for (comp in names(offsets)) {
    for (k in seq_len(n)) {
      out[[paste0(comp, k)]] <- as.vector(raw[offsets[[comp]] + k, , ])
    }
  }
  out[]
}

# Stochastic adjacency-list SIR (discrete-time tau-leaping) — the stochastic
# counterpart of running run_det_cd(..., sparse=TRUE). Takes a binary contact
# matrix, converts it to an adjacency list, then runs the stoch_mod_adj.R
# model via dust2 with n_particles = n_sim.
#
# Defaults follow the individual-level convention used by run_mean_field:
# one individual per node (N = rep(1, n)). Pass an explicit S_ini / I_ini for
# grouped use.
run_stoch_adj <- function(contact_matrix, beta, t, I_ini,
                          N=NULL, susceptibility=NULL, transmissibility=NULL,
                          gamma=1/3, waning=0, dt=0.1,
                          timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n   <- nrow(contact_matrix)
  adj <- contact_matrix_to_adj(contact_matrix)

  if (is.null(N))                N                <- rep(1L, n)
  if (is.null(susceptibility))   susceptibility   <- rep(1, n)
  if (is.null(transmissibility)) transmissibility <- rep(1, n)

  beta_scalar <- if (length(beta) == 1L) beta else mean(beta)

  params <- list(
    n               = as.integer(n),
    max_degree      = as.integer(adj$max_degree),
    beta            = beta_scalar,
    gamma           = gamma,
    waning          = waning,
    neighbors       = adj$neighbors,
    mask            = adj$mask,
    susceptibility  = susceptibility,
    transmisibility = transmissibility,
    S_ini           = N - I_ini,
    I_ini           = I_ini
  )

  sys <- dust2::dust_system_create(stoch_model_adj, params,
                                   n_particles = n_sim,
                                   n_threads   = cores,
                                   dt          = dt,
                                   time        = 0)
  dust2::dust_system_set_state_initial(sys)
  raw <- dust2::dust_system_simulate(sys, timepoints)
  if (length(dim(raw)) == 2L) {
    raw <- array(raw, dim = c(dim(raw)[1L], 1L, dim(raw)[2L]))
  }
  # State layout (declaration order in stoch_mod_adj.R):
  #   S[1..n], I[1..n], R[1..n], C[1..n]
  n_t <- length(timepoints)
  out <- data.table(
    time = rep(timepoints,     each  = n_sim),
    sim  = rep(seq_len(n_sim), times = n_t)
  )
  offsets <- c(S=0L, I=n, R=2L*n, C=3L*n)
  for (comp in names(offsets)) {
    for (k in seq_len(n)) {
      out[[paste0(comp, k)]] <- as.vector(raw[offsets[[comp]] + k, , ])
    }
  }
  out[]
}

# Individual-level stochastic SIR backed by stoch_ind.R. Takes a weighted
# contact matrix (entries are contact strengths; 0 means no contact). I_ini
# is a 0/1 vector marking initially infected individuals. Output columns per
# (sim, time): S1..n, I1..n, R1..n, foi1..n — foi[k] at recorded time t is
# the FOI that drove person k's transition over the previous step.
run_stoch_ind <- function(contact_matrix, beta, t, I_ini,
                          susceptibility=NULL, transmissibility=NULL,
                          gamma=1/3, dt=0.1,
                          timepoints=seq(0, t, 1),
                          n_sim=100, cores=10) {
  n   <- nrow(contact_matrix)
  adj <- contact_matrix_to_weighted_adj(contact_matrix)

  if (is.null(susceptibility))   susceptibility   <- rep(1, n)
  if (is.null(transmissibility)) transmissibility <- rep(1, n)

  I_ini_int <- as.integer(I_ini)
  S_ini_int <- 1L - I_ini_int

  params <- list(
    n               = as.integer(n),
    max_degree      = as.integer(adj$max_degree),
    beta            = beta,
    gamma           = gamma,
    neighbors       = adj$neighbors,
    weights         = adj$weights,
    susceptibility  = susceptibility,
    transmisibility = transmissibility,
    S_ini           = S_ini_int,
    I_ini           = I_ini_int
  )

  sys <- dust2::dust_system_create(stoch_model_ind, params,
                                   n_particles = n_sim,
                                   n_threads   = cores,
                                   dt          = dt,
                                   time        = 0)
  dust2::dust_system_set_state_initial(sys)
  raw <- dust2::dust_system_simulate(sys, timepoints)
  if (length(dim(raw)) == 2L) {
    raw <- array(raw, dim = c(dim(raw)[1L], 1L, dim(raw)[2L]))
  }
  # State layout (declaration order in stoch_ind.R):
  #   S[1..n], I[1..n], R[1..n], foi[1..n]
  n_t <- length(timepoints)
  out <- data.table(
    time = rep(timepoints,     each  = n_sim),
    sim  = rep(seq_len(n_sim), times = n_t)
  )
  offsets <- c(S=0L, I=n, R=2L*n, foi=3L*n)
  for (comp in names(offsets)) {
    for (k in seq_len(n)) {
      out[[paste0(comp, k)]] <- as.vector(raw[offsets[[comp]] + k, , ])
    }
  }
  out[]
}

# Stochastic counterpart of run_mean_field (det_model.R): same signature/
# defaults, backed by run_stoch_adj. Builds a Pareto contact network if c_ij
# is not given, picks a vaccination set if vac is not given, then runs the
# adjacency-list stochastic SIR with susceptibility = alpha for vaccinated
# nodes and 1 otherwise. Returns a list with $sum (per-sim/time vac/unvac/CRR)
# and $full (the underlying per-node data.table from run_stoch_adj).
run_stoch_network <- function(beta=1, N=100, pl_alpha=3,
                              susceptibility=c(1, 1), t=100,
                              vac_frac=0.5, vac=NULL, gamma=1/3,
                              c_ij=NULL, k_mean=6,
                              dt=0.1, timepoints=seq(0, t, 1),I_ini=2,
                              n_sim=100, cores=10) {
  if (is.null(c_ij))  c_ij <- get_conact_matrix_pl(N, pl_alpha, mean_k=k_mean)
  if (is.null(vac))   vac  <- sample(seq_len(N), vac_frac * N)

  # susceptibility = c(control, vaccinated)
  susept <- rep(susceptibility[1], N)
  susept[vac] <- susceptibility[2]
  non_vac <- setdiff(seq_len(N), vac)

  I_ini <- c(rep(1, I_ini), rep(0L, N - I_ini))

  full <- run_stoch_adj(c_ij, beta = N * beta / k_mean, t = t, I_ini = I_ini,
                        susceptibility = susept, gamma = gamma,
                        dt = dt, timepoints = timepoints,
                        n_sim = n_sim, cores = cores)

  vac_cols   <- paste0("C", vac)
  unvac_cols <- paste0("C", non_vac)
  sum_df <- full[, .(time, sim,
                     C2   = rowSums(.SD[, vac_cols,   with=FALSE]),
                     C1 = rowSums(.SD[, unvac_cols, with=FALSE]))]
#  sum_df[, CRR := (vac / (vac_frac * N)) / (unvac / ((1 - vac_frac) * N))]

  return(sum_df)#list(sum = sum_df, full = full)
}

# Generate directed contact events for all individual pairs.
# Rate (i->j): max_sus[k_j] * beta * mixing_matrix[k_j, k_i] / N_total.
# Returns a list of (from, to, t, unif) sorted by t.
generate_sir_contacts <- function(N_total, group_of, mixing_matrix, beta, max_sus, T) {
  contacts <- list()
  for (from in seq_len(N_total)) {
    k_f <- group_of[from]
    for (to in seq_len(N_total)) {
      if (from == to) next
      k_t  <- group_of[to]
      rate <- max_sus[k_t] * beta * mixing_matrix[k_t, k_f] / N_total
      if (rate <= 0) next
      t_curr <- 0
      while (TRUE) {
        t_curr <- t_curr + rexp(1, rate=rate)
        if (t_curr >= T) break
        contacts <- c(contacts, list(list(from=from, to=to, t=t_curr, unif=runif(1))))
      }
    }
  }
  contacts[order(sapply(contacts, `[[`, "t"))]
}

# Apply shared contacts to one scenario.
# p_sus[i]: per-individual acceptance probability = sus[k_i] / max_sus[k_i].
# recovery_times[i]: pre-generated recovery duration for individual i.
# Returns inf_times (Inf if never infected).
run_sir_events <- function(N_total, I_ini_expanded, p_sus, contacts, recovery_times) {
  inf_times <- rep(Inf, N_total)
  rec_times <- rep(Inf, N_total)
  for (i in which(as.logical(I_ini_expanded))) {
    inf_times[i] <- 0
    rec_times[i] <- recovery_times[i]
  }
  for (ev in contacts) {
    from <- ev$from; to <- ev$to; t_ev <- ev$t
    if (inf_times[from] > t_ev || rec_times[from] <= t_ev) next
    if (is.finite(inf_times[to])) next
    if (ev$unif < p_sus[to]) {
      inf_times[to] <- t_ev
      rec_times[to] <- t_ev + recovery_times[to]
    }
  }
  inf_times
}

# Coupled individual-level SIR: both scenarios share the same contact events and
# per-individual recovery durations; susceptibility differences are handled by thinning.
# sus_x / sus_z: per-group susceptibility (can be > 1). Returns Cx1..n, Cz1..n, time, sim.
run_coupled_sir <- function(beta, N, mixing_matrix, t, I_ini, sus_x, sus_z,
                             gamma=1/3, timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n_groups <- length(N)
  group_of <- rep(seq_len(n_groups), N)
  N_total  <- sum(N)
  max_sus  <- pmax(sus_x, sus_z)

  p_sus_x_ind <- (sus_x / max_sus)[group_of]
  p_sus_z_ind <- (sus_z / max_sus)[group_of]

  I_ini_expanded <- rep(0L, N_total)
  cumN <- cumsum(c(0L, N))
  for (k in seq_len(n_groups)) {
    if (I_ini[k] > 0L)
      I_ini_expanded[(cumN[k]+1L):min(cumN[k]+I_ini[k], cumN[k+1])] <- 1L
  }

  idx_of   <- lapply(seq_len(n_groups), function(k) which(group_of == k))
  Cx_names <- paste0("Cx", seq_len(n_groups))
  Cz_names <- paste0("Cz", seq_len(n_groups))

  run_one <- function(sim_id) {
    contacts       <- generate_sir_contacts(N_total, group_of, mixing_matrix, beta, max_sus, t)
    recovery_times <- rexp(N_total, rate=gamma)

    inf_x <- run_sir_events(N_total, I_ini_expanded, p_sus_x_ind, contacts, recovery_times)
    inf_z <- run_sir_events(N_total, I_ini_expanded, p_sus_z_ind, contacts, recovery_times)

    rows <- lapply(timepoints, function(tp) {
      row <- list(time=tp, sim=sim_id)
      for (k in seq_len(n_groups)) {
        row[[Cx_names[k]]] <- sum(inf_x[idx_of[[k]]] <= tp)
        row[[Cz_names[k]]] <- sum(inf_z[idx_of[[k]]] <= tp)
      }
      row
    })
    rbindlist(lapply(rows, as.data.frame))
  }

  rbindlist(parallel::mclapply(seq_len(n_sim), run_one, mc.cores=cores))
}

run_coupled_frailty_sir <- function(alpha, sd, beta=1, R=NULL, f=0.5, N=200, t=30,
                                     n_frailty=10, gamma=1/2, vac_counts_x=NULL, vac_counts_z=NULL,
                                     I_ini_total=1, timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  if (!is.null(R)) beta <- get_beta(R, alpha, sd, f=f, N=N, n_frailty=n_frailty, gamma=gamma)

  fr      <- get_frailty(sd=sd, n=n_frailty)
  frailty <- exp(2.5 * fr$x)
  n_total <- round(2 * N * fr$p)
  if (is.null(vac_counts_x)) vac_counts_x <- round(f * n_total)
  if (is.null(vac_counts_z)) vac_counts_z <- round(f * n_total)
  vac_counts_x <- expand_vac_counts(vac_counts_x, n_total)
  vac_counts_z <- expand_vac_counts(vac_counts_z, n_total)

  N_total <- sum(n_total)
  cumN    <- cumsum(c(0L, n_total))

  sus_x   <- numeric(N_total); sus_z <- numeric(N_total)
  vac_x   <- logical(N_total); vac_z <- logical(N_total)
  for (k in seq_len(n_frailty)) {
    idx  <- (cumN[k]+1L):cumN[k+1L]
    n_k  <- n_total[k]
    is_x <- seq_len(n_k) <= vac_counts_x[k]
    is_z <- seq_len(n_k) <= vac_counts_z[k]
    sus_x[idx] <- ifelse(is_x, alpha * frailty[k], frailty[k])
    sus_z[idx] <- ifelse(is_z, alpha * frailty[k], frailty[k])
    vac_x[idx] <- is_x; vac_z[idx] <- is_z
  }

  N_vac_x   <- sum(vac_counts_x); N_unvac_x <- N_total - N_vac_x
  N_vac_z   <- sum(vac_counts_z); N_unvac_z <- N_total - N_vac_z

  max_sus     <- pmax(sus_x, sus_z)
  p_sus_x_ind <- sus_x / max_sus
  p_sus_z_ind <- sus_z / max_sus

  # Seed: first unvaccinated individual in the largest frailty group (shared between scenarios)
  I_ini_expanded <- rep(0L, N_total)
  k_seed <- which.max(n_total)
  seed_pool <- which(!vac_x & seq_len(N_total) %in% ((cumN[k_seed]+1L):cumN[k_seed+1L]))
  if (length(seed_pool) > 0)
    I_ini_expanded[seed_pool[seq_len(min(I_ini_total, length(seed_pool)))]] <- 1L

  # Individual-level contacts with homogeneous mixing consistent with the frailty group model.
  # The frailty group model uses mm = 1/n_groups; replicating at individual level requires
  # mm_ind[i,j] = 1/n_groups so that generate_sir_contacts gives rate = max_sus * beta / (n_groups * N_total).
  group_of <- seq_len(N_total)
  mm_ind   <- matrix(1 / (2 * n_frailty), nrow=N_total, ncol=N_total)

  run_one <- function(sim_id) {
    contacts       <- generate_sir_contacts(N_total, group_of, mm_ind, beta, max_sus, t)
    recovery_times <- rexp(N_total, rate=gamma)

    inf_x <- run_sir_events(N_total, I_ini_expanded, p_sus_x_ind, contacts, recovery_times)
    inf_z <- run_sir_events(N_total, I_ini_expanded, p_sus_z_ind, contacts, recovery_times)

    rbindlist(lapply(timepoints, function(tp) {
      data.frame(time=tp, sim=sim_id,
                 vac_x   = sum(inf_x[vac_x]  <= tp),
                 unvac_x = sum(inf_x[!vac_x] <= tp),
                 CRR_x   = (sum(inf_x[vac_x] <= tp) / N_vac_x) /
                            (sum(inf_x[!vac_x] <= tp) / N_unvac_x),
                 vac_z   = sum(inf_z[vac_z]  <= tp),
                 unvac_z = sum(inf_z[!vac_z] <= tp),
                 CRR_z   = (sum(inf_z[vac_z] <= tp) / N_vac_z) /
                            (sum(inf_z[!vac_z] <= tp) / N_unvac_z))
    }))
  }

  rbindlist(parallel::mclapply(seq_len(n_sim), run_one, mc.cores=cores))
}

run_stoch_linear <- function(beta, N, t, susceptibility=NULL,
                             timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n <- length(N)
  if (is.null(susceptibility)) susceptibility <- rep(1, n)

  S_names <- paste0("S", 1:n)
  C_names <- paste0("C", 1:n)

  transitions <- lapply(1:n, function(i)
    setNames(c(-1L, 1L), c(S_names[i], C_names[i]))
  )

  rate_func <- local({
    sus <- susceptibility
    function(x, p, t) pmax(0, sus * beta * x[S_names])
  })

  init <- setNames(c(rbind(N, rep(0, n))), c(rbind(S_names, C_names)))

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions, rate_func, list(), tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

# odin2/dust2 equivalent of run_stoch_linear. Same R-level signature; uses
# stoch_linear.R via dust2 with n_particles = n_sim. See the note on `dt` for
# run_stoch_cd_dust — shrink dt when susceptibility * beta is large.
run_stoch_linear_dust <- function(beta, N, t, susceptibility=NULL,
                                  dt=0.1, timepoints=seq(0, t, 1),
                                  n_sim=100, cores=10) {
  n <- length(N)
  if (is.null(susceptibility)) susceptibility <- rep(1, n)

  params <- list(
    n              = as.integer(n),
    beta           = beta,
    S_ini          = N,
    susceptibility = susceptibility
  )

  sys <- dust2::dust_system_create(stoch_model_linear, params,
                                   n_particles = n_sim,
                                   n_threads   = cores,
                                   dt          = dt,
                                   time        = 0)
  dust2::dust_system_set_state_initial(sys)
  raw <- dust2::dust_system_simulate(sys, timepoints)
  if (length(dim(raw)) == 2L) {
    raw <- array(raw, dim = c(dim(raw)[1L], 1L, dim(raw)[2L]))
  }
  # State layout (declaration order in stoch_linear.R): S[1..n], C[1..n]
  n_t <- length(timepoints)
  out <- data.table(
    time = rep(timepoints,     each  = n_sim),
    sim  = rep(seq_len(n_sim), times = n_t)
  )
  offsets <- c(S=0L, C=n)
  for (comp in names(offsets)) {
    for (k in seq_len(n)) {
      out[[paste0(comp, k)]] <- as.vector(raw[offsets[[comp]] + k, , ])
    }
  }
  out[]
}

run_stoch_frailty_linear <- function(alpha, sd, beta=1, f=0.5, N=1000, t=100,
                                      n_frailty=100, vac_counts=NULL,
                                      timepoints=seq(0, t, 1), n_sim=100, cores=10,
                                      method=c("ctmc", "dust"), dt=0.1) {
  method <- match.arg(method)
  fr      <- get_frailty(sd=sd, n=n_frailty)
  frailty <- exp(2.5 * fr$x)
  n_total <- round(2 * N * fr$p)
  if (is.null(vac_counts)) vac_counts <- round(f * n_total)

  n_groups       <- 2 * n_frailty
  N_groups       <- c(n_total - vac_counts, vac_counts)
  susceptibility <- c(frailty, alpha * frailty)

  raw <- if (method == "ctmc") {
    run_stoch_linear(beta=beta, N=N_groups, t=t, susceptibility=susceptibility,
                     timepoints=timepoints, n_sim=n_sim, cores=cores)
  } else {
    run_stoch_linear_dust(beta=beta, N=N_groups, t=t, susceptibility=susceptibility,
                          dt=dt, timepoints=timepoints, n_sim=n_sim, cores=cores)
  }

  N_vac   <- sum(vac_counts)
  N_unvac <- sum(n_total - vac_counts)
  raw[, vac   := rowSums(.SD), .SDcols=paste0("C", (n_frailty+1):n_groups)]
  raw[, unvac := rowSums(.SD), .SDcols=paste0("C", 1:n_frailty)]
  raw[, CRR   := (vac / N_vac) / (unvac / N_unvac)]
  raw[, .(time, sim, vac, unvac, CRR)]
}

# Generate Poisson exposure events for N individuals at per-individual rates.
# Returns data.frame(who, t, unif) sorted by t.
generate_linear_events <- function(N, rates, T) {
  events <- list()
  for (i in seq_len(N)) {
    if (rates[i] <= 0) next
    t_curr <- 0
    while (TRUE) {
      t_curr <- t_curr + rexp(1, rate=rates[i])
      if (t_curr >= T) break
      events <- c(events, list(c(who=i, t=t_curr, unif=runif(1))))
    }
  }
  if (length(events) == 0)
    return(data.frame(who=integer(), t=numeric(), unif=numeric()))
  df <- as.data.frame(do.call(rbind, events))
  df[order(df$t), ]
}

# Coupled linear model: two scenarios (X and Z) share the same Poisson event stream.
# Events are generated at rate max(sus_x[k], sus_z[k]) * beta per individual;
# each event is accepted for scenario X with probability sus_x[k] / max_sus[k],
# and for Z with probability sus_z[k] / max_sus[k].
# sus_x / sus_z: per-group susceptibility (can be > 1).
# Returns columns Cx1..n, Cz1..n, time, sim.
run_coupled_linear <- function(beta, N, t, sus_x, sus_z,
                                timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n        <- length(N)
  max_sus  <- pmax(sus_x, sus_z)
  group_of <- rep(seq_len(n), N)
  N_total  <- sum(N)

  # Per-individual Poisson rates and acceptance probabilities
  rates <- max_sus[group_of] * beta
  p_x   <- (sus_x / max_sus)[group_of]
  p_z   <- (sus_z / max_sus)[group_of]

  # Precompute per-group individual indices
  idx_of <- lapply(seq_len(n), function(k) which(group_of == k))

  Cx_names <- paste0("Cx", seq_len(n))
  Cz_names <- paste0("Cz", seq_len(n))

  run_one <- function(sim_id) {
    events <- generate_linear_events(N_total, rates, t)
    inf_x  <- rep(Inf, N_total)
    inf_z  <- rep(Inf, N_total)

    for (j in seq_len(nrow(events))) {
      i <- events$who[j]
      u <- events$unif[j]
      if (!is.finite(inf_x[i]) && u < p_x[i]) inf_x[i] <- events$t[j]
      if (!is.finite(inf_z[i]) && u < p_z[i]) inf_z[i] <- events$t[j]
    }

    rows <- lapply(timepoints, function(tp) {
      row <- list(time=tp, sim=sim_id)
      for (k in seq_len(n)) {
        row[[Cx_names[k]]] <- sum(inf_x[idx_of[[k]]] <= tp)
        row[[Cz_names[k]]] <- sum(inf_z[idx_of[[k]]] <= tp)
      }
      row
    })
    rbindlist(lapply(rows, as.data.frame))
  }

  rbindlist(parallel::mclapply(seq_len(n_sim), run_one, mc.cores=cores))
}

# If vac_counts is a scalar, distribute it proportionally across frailty groups by size.
expand_vac_counts <- function(vac_counts, n_total) {
  if (length(vac_counts) == 1L)
    pmin(round(vac_counts * n_total / sum(n_total)), n_total)
  else
    vac_counts
}

run_coupled_frailty_linear <- function(alpha, sd, beta=1, f=0.5, N=1000, t=100,
                                        n_frailty=100, vac_counts_x=NULL, vac_counts_z=NULL,
                                        timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  fr      <- get_frailty(sd=sd, n=n_frailty)
  frailty <- exp(2.5 * fr$x)
  n_total <- round(2 * N * fr$p)
  if (is.null(vac_counts_x)) vac_counts_x <- round(f * n_total)
  if (is.null(vac_counts_z)) vac_counts_z <- round(f * n_total)
  vac_counts_x <- expand_vac_counts(vac_counts_x, n_total)
  vac_counts_z <- expand_vac_counts(vac_counts_z, n_total)

  N_total <- sum(n_total)
  cumN    <- cumsum(c(0L, n_total))

  sus_x    <- numeric(N_total); sus_z    <- numeric(N_total)
  is_vac_x <- logical(N_total); is_vac_z <- logical(N_total)
  for (k in seq_len(n_frailty)) {
    idx  <- (cumN[k]+1L):cumN[k+1L]
    n_k  <- n_total[k]
    ix   <- seq_len(n_k) <= vac_counts_x[k]
    iz   <- seq_len(n_k) <= vac_counts_z[k]
    sus_x[idx]    <- ifelse(ix, alpha * frailty[k], frailty[k])
    sus_z[idx]    <- ifelse(iz, alpha * frailty[k], frailty[k])
    is_vac_x[idx] <- ix; is_vac_z[idx] <- iz
  }

  N_vac_x <- sum(vac_counts_x); N_unvac_x <- N_total - N_vac_x
  N_vac_z <- sum(vac_counts_z); N_unvac_z <- N_total - N_vac_z

  raw <- run_coupled_linear(beta=beta, N=rep(1L, N_total), t=t,
                             sus_x=sus_x, sus_z=sus_z,
                             timepoints=timepoints, n_sim=n_sim, cores=cores)

  Cx_cols <- paste0("Cx", seq_len(N_total))
  Cz_cols <- paste0("Cz", seq_len(N_total))
  agg <- function(cols) if (length(cols)) rowSums(raw[, cols, with=FALSE]) else 0L
  raw[, vac_x   := agg(Cx_cols[ is_vac_x])]
  raw[, unvac_x := agg(Cx_cols[!is_vac_x])]
  raw[, vac_z   := agg(Cz_cols[ is_vac_z])]
  raw[, unvac_z := agg(Cz_cols[!is_vac_z])]
  raw[, CRR_x   := (vac_x / N_vac_x) / (unvac_x / N_unvac_x)]
  raw[, CRR_z   := (vac_z / N_vac_z) / (unvac_z / N_unvac_z)]
  raw[, .(time, sim, vac_x, unvac_x, CRR_x, vac_z, unvac_z, CRR_z)]
}

# `sd_trans` adds a per-group transmissibility frailty that is perfectly rank-
# correlated with the susceptibility frailty: each frailty bin k has a
# susceptibility Beta-rank, and trans uses the matching quantile mapped
# through Beta(sd_trans). With sd_trans = sd this reproduces the
# susceptibility values exactly; with sd_trans = 0 there is no trans
# heterogeneity. Vaccination is taken to act on susceptibility only (alpha).
# get_beta sees sd_trans and incorporates it into the NGM, so R calibration
# is consistent with the simulator.
run_stoch_frailty_cd <- function(sd, sd_trans=0, beta=1, R=NULL, f=0.5, N=1000, t=100,
                                  n_frailty=100, gamma=1/2, vac_counts=NULL,
                                  I_ini_total=1, timepoints=seq(0, t, 1), n_sim=100, cores=10,
                                  method=c("ctmc", "dust"), dt=0.1, susceptibility=c(1,0.5)) {
  method <- match.arg(method)
  alpha <- susceptibility[2] 
  if (!is.null(R)) beta <- get_beta(R, alpha, sd, sd_trans=sd_trans, f=f, N=N, n_frailty=n_frailty, gamma=gamma)

  # Population-grid frailty: use whichever of sd/sd_trans is > 0 for the bin
  # weighting; if both > 0 we key on `sd` (consistent with the prior behaviour
  # when sd_trans was not yet a knob).
  sd_pop <- if (sd > 0) sd else sd_trans
  fr     <- get_frailty(sd=sd_pop, n=n_frailty)
  n_total <- round(2 * N * fr$p)
  if (is.null(vac_counts)) vac_counts <- round(f * n_total)

  # Per-bin susceptibility frailty. When sd = 0 there is no sus heterogeneity.
  frailty <- if (sd > 0) exp(2.5 * fr$x) else rep(1, n_frailty)

  # Per-bin transmissibility frailty, rank-correlated with sus. If sd_trans
  # matches the population sd_pop the bins already sit at the right values;
  # otherwise map ranks through the matching Beta inverse CDF.
  trans_frailty <- if (sd_trans > 0) {
    if (sd_trans == sd_pop) {
      exp(2.5 * fr$x)
    } else {
      cf_pop        <- (0.25 / sd_pop^2)   - 1
      cf_t          <- (0.25 / sd_trans^2) - 1
      ranks         <- pbeta(fr$x, 0.5 * cf_pop, 0.5 * cf_pop)
      exp(2.5 * qbeta(ranks, 0.5 * cf_t, 0.5 * cf_t))
    }
  } else {
    rep(1, n_frailty)
  }

  n_groups       <- 2 * n_frailty
  N_groups       <- c(n_total - vac_counts, vac_counts)
  susceptibility <- c(frailty, alpha * frailty)
  transmissibility <- c(trans_frailty, trans_frailty)  # vaccine acts on sus only
  mm             <- matrix(1, nrow=n_groups, ncol=n_groups) / n_groups
  
  # Spread I_ini_total seeds proportionally across vac and unvac populations,
  # then within each across frailty bins. Two-step rounding keeps the
  # aggregate vac/unvac split true to f (e.g. f=0.5 gives a 50/50 seed split
  # by population, not by which bin happened to win the floor rounding).
  if (I_ini_total > sum(N_groups))
    stop(sprintf("I_ini_total (%d) exceeds total population (%d).",
                 I_ini_total, sum(N_groups)))
  spread <- function(total, sizes) {
    if (total == 0L) return(integer(length(sizes)))
    target <- total * sizes / sum(sizes)
    ini    <- pmin(floor(target), sizes)
    rem    <- total - sum(ini)
    if (rem > 0L) {
      for (idx in order(target - ini, decreasing = TRUE)) {
        if (rem == 0L) break
        if (ini[idx] < sizes[idx]) {
          ini[idx] <- ini[idx] + 1L
          rem      <- rem - 1L
        }
      }
    }
    as.integer(ini)
  }
  N_unvac     <- N_groups[seq_len(n_frailty)]
  N_vac       <- N_groups[(n_frailty + 1L):n_groups]
  total_unvac <- sum(N_unvac); total_vac <- sum(N_vac)
  ut          <- I_ini_total * total_unvac / (total_unvac + total_vac)
  unvac_seeds <- min(floor(ut), total_unvac)
  vac_seeds   <- I_ini_total - unvac_seeds
  # If the rounding "owes" the unvac side the half-seed back, take it.
  if (vac_seeds > total_vac) {
    deficit     <- vac_seeds - total_vac
    vac_seeds   <- as.integer(total_vac)
    unvac_seeds <- unvac_seeds + deficit
  } else if (ut - unvac_seeds > 0.5 && unvac_seeds < total_unvac && vac_seeds > 0L) {
    unvac_seeds <- unvac_seeds + 1L; vac_seeds <- vac_seeds - 1L
  }
  I_ini <- c(spread(as.integer(unvac_seeds), N_unvac),
             spread(as.integer(vac_seeds),   N_vac))

  raw <- if (method == "ctmc") {
    run_stoch_cd_ctmc(mm, beta=beta, N=N_groups, t=t, I_ini=I_ini,
                      susceptibility=susceptibility,
                      transmissibility=transmissibility, gamma=gamma,
                      timepoints=timepoints, n_sim=n_sim, cores=cores)
  } else {
    run_stoch_cd_dust(mm, beta=beta, N=N_groups, t=t, I_ini=I_ini,
                      susceptibility=susceptibility,
                      transmissibility=transmissibility, gamma=gamma,
                      dt=dt, timepoints=timepoints, n_sim=n_sim, cores=cores)
  }

  N_vac   <- sum(vac_counts)
  N_unvac <- sum(n_total - vac_counts)
  raw[, vac   := rowSums(.SD), .SDcols=paste0("C", (n_frailty+1):n_groups)]
  raw[, unvac := rowSums(.SD), .SDcols=paste0("C", 1:n_frailty)]
  raw[, CRR   := (vac / N_vac) / (unvac / N_unvac)]
  raw[, C1:=unvac]
  raw[, C2:=vac]
  raw[, .(time, sim, C1, C2, vac, unvac, CRR)]
}

regularise <- function(df, timepoints) {
  df_reg <- dplyr::bind_rows(df, data.frame(time=timepoints)) |> dplyr::arrange(time)
  for (col in setdiff(colnames(df_reg), c("time", "sim"))) {
    df_reg[[col]] <- zoo::na.locf(df_reg[[col]], na.rm=FALSE)
  }
  dplyr::distinct(dplyr::filter(df_reg, time %in% timepoints))
}



get_lik <- function(par, mod, X_vac, X_cont, delta= 3, eta=1e-10){

  if(any(par <= 0)) return(-Inf)
  #print(exp(par))
  a <- mod(beta=par[1], susceptibility=c(1, par[2]))
  max_t <- max(a$time)
  setDT(a)
  lik <- nrow(a[time==max_t & abs(C1 - X_cont) < delta& abs(C2 - X_vac) < delta,])/nrow(a[time==max_t,])
  log_lik <- log((1-eta)*lik + eta/(nrow(a[time==max_t,])))
  log_prio <- ifelse(all(par>0 & par<3), 0, -Inf)
 # print(log_lik)
  return(log_lik + log_prio)
}

# Chunked Metropolis: run `n` steps in `progress_every`-sized batches,
# concatenating samples and printing a one-line summary per chunk so long
# fits aren't silent. Each chunk continues from the previous final state,
# so the resulting chain is identical to a single nbatch=n call modulo RNG
# stream resets between chunks.
run_metrop_chunked <- function(lik, initial, n, scale, progress_every=100, label="") {
  batches      <- list()
  accept_sum   <- 0
  accept_count <- 0
  current      <- as.numeric(initial)
  done         <- 0
  while (done < n) {
    chunk        <- min(progress_every, n - done)
    res          <- mcmc::metrop(lik, initial=current, nbatch=chunk, scale=scale)
    batches[[length(batches)+1]] <- res$batch
    accept_sum   <- accept_sum   + res$accept * chunk
    accept_count <- accept_count + chunk
    current      <- as.numeric(tail(res$batch, 1))
    done         <- done + chunk
    message(sprintf("[%s] %d/%d  accept=%.3f  current=(%s)",
                    label, done, n, accept_sum/accept_count,
                    paste(sprintf("%.4f", current), collapse=", ")))
  }
  list(batch  = do.call(rbind, batches),
       accept = accept_sum / accept_count,
       final  = current)
}

fit_mod <- function(mod, X_cont=NULL, X_vac=NULL, beta_ini=1, alpha_ini=0.5, n=500, burn_in=100, scale=matrix(c(0.02, 0.02, 0.05, 0.05), nrow=2), progress_every=100) {
   lik <- purrr::partial(get_lik, mod=mod, X_vac=X_vac, X_cont=X_cont)

   res_first <- run_metrop_chunked(lik, initial=c(beta=beta_ini, alpha=alpha_ini),
                                   n=burn_in, scale=scale,
                                   progress_every=progress_every, label="burn-in")
    Sigma_post <- cov(res_first$batch)
    current <- tail(res_first$batch, 1)
    d <- length(current)
    new_cov <- t(chol(2.38^2/d * Sigma_post))
    print("Covariance matrix:")
    print(new_cov)
   res <- run_metrop_chunked(lik, initial=current, n=n, scale=new_cov,
                             progress_every=progress_every, label="sample")

   print(" ### SUMMARY ### ")
   print("Acceptance rate:")
   print(res$accept)
   print("Effective sample size:")
   print(coda::effectiveSize(res$batch))

  return(res$batch)

}

get_lik_norm <- function(par, mod, X_vac, X_cont, cor_matrix){

  if(any(par <= 0)) return(-Inf)
  #print(exp(par))
  a <- mod(beta_day=par[1], susceptibility=c(1, par[2]))$main
  setDT(a)
  mu1 <- as.numeric(a[t==max(a$t),"unvac"])
  mu2 <- as.numeric(a[t==max(a$t), "vac"])
#  print(mu1)
 # print(mu2)
  #print(cor_matrix)
  #print(c(X_cont, X_vac))
  lik <- mvtnorm::dmvnorm(c(X_cont, X_vac), mean=c(mu1, mu2), sigma=cor_matrix)
  log_lik <- log(lik)
  log_prio <- ifelse(all(par>0 & par<3), 0, -Inf)
 # print(log_lik)
  return(log_lik + log_prio)
}

fit_mod_norm <- function(mod, X_cont=NULL, X_vac=NULL, beta_ini=1, alpha_ini=0.5,cor_matrix=matrix(c(10,10,10,10), nrow=2), n=500, burn_in=100, scale=matrix(c(0.02, 0.02, 0.05, 0.05), nrow=2), progress_every=100) {
   lik <- purrr::partial(get_lik_norm, mod=mod, X_vac=X_vac, X_cont=X_cont, cor_matrix=cor_matrix)
   res_first <- run_metrop_chunked(lik, initial=c(beta=beta_ini, alpha=alpha_ini),
                                   n=burn_in, scale=scale,
                                   progress_every=progress_every, label="burn-in")
    Sigma_post <- cov(res_first$batch)
    current <- tail(res_first$batch, 1)
    d <- length(current)
    new_cov <- t(chol(2.38^2/d * Sigma_post))
    print("Covariance matrix:")
    print(new_cov)
   res <- run_metrop_chunked(lik, initial=current, n=n, scale=new_cov,
                             progress_every=progress_every, label="sample")

   print(" ### SUMMARY ### ")
   print("Acceptance rate:")
   print(res$accept)
   print("Effective sample size:")
   print(coda::effectiveSize(res$batch))

  return(res$batch)
  
}
