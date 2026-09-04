library(adaptivetau)
library(data.table)

source("odin_cache.R")

# Event-driven network simulator (adj_to_csr, run_stoch_network_events,
# net_sir_compile). Sourcing only defines functions -- the C++ translation
# unit is compiled lazily on first use -- so this is free for callers that
# stay on the dust engine. Guarded so a machine without cpp11 can still load
# stoch_model.R and run everything except the "events" engine.
if (requireNamespace("cpp11", quietly = TRUE)) {
  local({
    f <- file.path(odin_project_root(), "net_sir_events.R")
    if (file.exists(f)) source(f)
  })
}

stoch_model_cd         <- odin_cached("stoch_sir.R")
stoch_model_adj        <- odin_cached("stoch_mod_adj.R")
stoch_model_adj_sparse <- odin_cached("stoch_mod_adj_sparse.R")
stoch_model_adj_vacfrac<- odin_cached("stoch_mod_adj_vacfrac.R")
stoch_model_linear     <- odin_cached("stoch_linear.R")
stoch_model_ind        <- odin_cached("stoch_ind.R")

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
                              timepoints=seq(0, t, 1), n_sim=100, cores=10,
                              seed=NULL) {
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
                                   time        = 0,
                                   seed        = seed)
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
                          timepoints=seq(0, t, 1), n_sim=100, cores=10,
                          seed=NULL, adj=NULL) {
  # `adj` lets the caller pass a PRE-BUILT adjacency list. contact_matrix_to_adj
  # is an R-level loop over every row, i.e. O(n^2) and slow (1.8 s at n = 5000),
  # and the contact matrix is fixed for a whole config -- so rebuilding it on
  # every simulator call dominated the network runtime. Build it once
  # (materialise_cfg) and pass it in.
  adj <- if (is.null(adj)) contact_matrix_to_adj(contact_matrix) else adj
  n   <- if (is.null(contact_matrix)) ncol(adj$neighbors) else nrow(contact_matrix)

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
                                   time        = 0,
                                   seed        = seed)
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


# Group-aggregate-output variant of run_stoch_adj. Same dynamics as the
# dense version but backed by stoch_mod_adj_sparse.R; only n_groups scalar
# S_g / I_g / R_g / C_g summaries are returned per (time, sim), so memory
# scales with n_groups instead of N. Use when N is large enough that the
# per-node output would be prohibitive.
#
#   adj                 : list(neighbors, mask, max_degree) from either
#                         contact_matrix_to_adj() or sample_pareto_adj().
#                         The full contact matrix is never needed.
#   groups              : named list of integer node-index vectors, one
#                         per group of interest (e.g. list(vac = vac_idx,
#                         control = ctrl_idx)). Groups can overlap.
#
# Returns a data.table with columns
#     time, sim, S_<group>, I_<group>, R_<group>, C_<group>
# for every named group.
run_stoch_adj_sparse <- function(adj, beta, t, I_ini, groups,
                                 N = NULL, susceptibility = NULL,
                                 transmissibility = NULL,
                                 gamma = 1 / 3, dt = 0.1,
                                 timepoints = seq(0, t, 1),
                                 n_sim = 100, cores = 10, seed = NULL) {
  n <- ncol(adj$neighbors)
  if (is.null(N))                N                <- rep(1L, n)
  if (is.null(susceptibility))   susceptibility   <- rep(1, n)
  if (is.null(transmissibility)) transmissibility <- rep(1, n)

  group_names <- names(groups)
  if (is.null(group_names) || any(!nzchar(group_names)))
    stop("`groups` must be a named list.")
  n_groups <- length(groups)
  G <- matrix(0L, nrow = n_groups, ncol = n)
  for (g in seq_len(n_groups)) {
    idx <- groups[[g]]
    if (any(idx < 1L | idx > n))
      stop("Node indices in group ", group_names[g], " are out of range.")
    G[g, idx] <- 1L
  }

  beta_scalar <- if (length(beta) == 1L) beta else mean(beta)

  params <- list(
    n               = as.integer(n),
    n_groups        = as.integer(n_groups),
    max_degree      = as.integer(adj$max_degree),
    beta            = beta_scalar,
    gamma           = gamma,
    neighbors       = adj$neighbors,
    mask            = adj$mask,
    susceptibility  = susceptibility,
    transmisibility = transmissibility,
    S_ini           = N - I_ini,
    I_ini           = I_ini,
    G               = G
  )

  sys <- dust2::dust_system_create(stoch_model_adj_sparse, params,
                                   n_particles = n_sim,
                                   n_threads   = cores,
                                   dt          = dt,
                                   time        = 0,
                                   seed        = seed)
  dust2::dust_system_set_state_initial(sys)

  # State layout (declaration order): S[1..n], I[1..n], R[1..n], C[1..n],
  # cumFOI[1..n], S_g[1..n_groups], I_g[..], R_g[..], C_g[..],
  # cumFOI_sum[1..n_groups], cumFOI_sumsq[1..n_groups]. Ask dust2 to
  # only return the aggregate block so the returned array stays
  # O(n_groups) per particle per timepoint.
  agg_start   <- 5L * n + 1L
  agg_indices <- seq(agg_start, agg_start + 6L * n_groups - 1L)
  raw <- dust2::dust_system_simulate(sys, timepoints, index_state = agg_indices)
  if (length(dim(raw)) == 2L)
    raw <- array(raw, dim = c(dim(raw)[1L], 1L, dim(raw)[2L]))

  n_t <- length(timepoints)
  out <- data.table(
    time = rep(timepoints,     each  = n_sim),
    sim  = rep(seq_len(n_sim), times = n_t)
  )
  offsets <- c(S = 0L, I = n_groups, R = 2L * n_groups, C = 3L * n_groups,
               cumFOI_sum = 4L * n_groups, cumFOI_sumsq = 5L * n_groups)
  for (comp in names(offsets)) {
    for (g in seq_len(n_groups)) {
      col <- paste0(comp, "_", group_names[g])
      out[[col]] <- as.vector(raw[offsets[[comp]] + g, , ])
    }
  }
  # Attach the per-group population sizes so callers can compute
  # per-group means and CV^2 directly (mean = sum / N_group; CV^2 =
  # sumsq * N_group / sum^2 - 1).
  attr(out, "N_group") <- setNames(lengths(groups), group_names)
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

# Which simulator run_stoch_network() uses under the hood:
#
#   "events" — exact event-driven C++ first-passage simulator
#              (net_sir_events.cpp). No dt discretisation bias, O(E log V)
#              per realisation, ~50-250x faster than "dust" at the dt the
#              pipeline needs. This is the default.
#   "dust"   — the original discrete-time tau-leaping model (stoch_mod_adj.R)
#              via run_stoch_adj(). Kept as the fallback: it is the only
#              engine that produces per-node S/I/R trajectories, and it is
#              what every result before this switch was produced with.
#
# Override globally with the environment variable (so a slurm array can be
# flipped without touching code) or the option:
#
#   EATE_NETWORK_ENGINE=dust Rscript run_fit_array.R
#   options(eate.network_engine = "dust")
#
# The env var wins; per-call `engine=` wins over both.
network_engine <- function() {
  e <- Sys.getenv("EATE_NETWORK_ENGINE", "")
  if (!nzchar(e)) e <- getOption("eate.network_engine", "events")
  match.arg(tolower(e), c("events", "dust"))
}

# Compile the event-driven simulator. It is a single translation unit and the
# build is cached for the session; call this eagerly in the parent process
# before any fork, so the mclapply workers inherit a loaded DLL instead of
# each triggering a build of their own.
ensure_net_sir_events <- function() {
  if (!exists("run_stoch_network_events", mode = "function"))
    stop("the \"events\" network engine needs net_sir_events.R and the 'cpp11' ",
         "package; install cpp11 or select the dust engine with ",
         "EATE_NETWORK_ENGINE=dust")
  net_sir_compile()
  invisible(TRUE)
}

# Stochastic counterpart of run_mean_field (det_model.R): same signature/
# defaults. Builds a Pareto contact network if c_ij is not given, picks a
# vaccination set if vac is not given, then runs the network SIR with
# susceptibility = alpha for vaccinated nodes and 1 otherwise. Returns a
# per-(sim, time) data.table of C1 (unvaccinated cases) and C2 (vaccinated).
#
# `engine` selects the simulator — see network_engine() above. Both engines
# simulate the same continuous-time process and return the same columns on
# the same scale; "dust" additionally honours `dt`, which "events" ignores
# because it has no time discretisation. `csr` is the event engine's
# pre-built graph (adj_to_csr); pass it for the same reason `adj` is passed
# to run_stoch_adj — so a fit does not rebuild it on every call.
run_stoch_network <- function(beta=1, N=100, pl_alpha=3,
                              susceptibility=c(1, 1), t=100,
                              vac_frac=0.5, vac=NULL, gamma=1/3,
                              c_ij=NULL, k_mean=6,
                              dt=0.1, timepoints=seq(0, t, 1),I_ini=2,
                              n_sim=100, cores=10, seed=NULL, adj=NULL,
                              engine=network_engine(), csr=NULL) {
  engine <- match.arg(engine, c("events", "dust"))
  if (is.null(c_ij) && is.null(adj) && is.null(csr))
    c_ij <- get_conact_matrix_pl(N, pl_alpha, mean_k=k_mean)
  if (is.null(vac))   vac  <- sample(seq_len(N), vac_frac * N)

  if (engine == "events") {
    ensure_net_sir_events()
    # beta is scaled inside run_stoch_network_events (same N/k_mean factor as
    # the dust branch below), so it takes the user-scale beta unchanged.
    return(run_stoch_network_events(
      beta = beta, N = N, susceptibility = susceptibility, t = t,
      vac = vac, csr = csr, adj = adj, c_ij = c_ij, gamma = gamma,
      timepoints = timepoints, I_ini = I_ini, n_sim = n_sim,
      seed = seed, k_mean = k_mean, cores = cores))
  }

  # susceptibility = c(control, vaccinated)
  susept <- rep(susceptibility[1], N)
  susept[vac] <- susceptibility[2]
  non_vac <- setdiff(seq_len(N), vac)

  I_ini <- c(rep(1, I_ini), rep(0L, N - I_ini))

  full <- run_stoch_adj(c_ij, beta = N * beta / k_mean, t = t, I_ini = I_ini,
                        susceptibility = susept, gamma = gamma,
                        dt = dt, timepoints = timepoints,
                        n_sim = n_sim, cores = cores, seed = seed, adj = adj)

  vac_cols   <- paste0("C", vac)
  unvac_cols <- paste0("C", non_vac)
  # C1 before C2, matching the "events" engine, so the two are comparable
  # column-for-column and not just by name.
  sum_df <- full[, .(time, sim,
                     C1 = rowSums(.SD[, unvac_cols, with=FALSE]),
                     C2 = rowSums(.SD[, vac_cols,   with=FALSE]))]
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
  frailty <- { raw <- exp(2.5 * fr$x); raw / sum(fr$p * raw) }
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
                                  n_sim=100, cores=10, seed=NULL) {
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
                                   time        = 0,
                                   seed        = seed)
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
  frailty <- { raw <- exp(2.5 * fr$x); raw / sum(fr$p * raw) }
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
  frailty <- { raw <- exp(2.5 * fr$x); raw / sum(fr$p * raw) }
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
                                  method=c("ctmc", "dust"), dt=0.1, susceptibility=c(1,0.5),
                                  frailty_amp=2.5, seed=NULL) {
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
  # Both arrays are normalised to population-mean 1 (sum(fr$p * x) = 1) so
  # heterogeneity contributes only via the NGM variance, not via a shift in
  # the mean transmission rate — beta therefore retains the R0 = beta/gamma
  # interpretation in the homogeneous limit.
  frailty <- if (sd > 0) {
    raw <- exp(frailty_amp * fr$x); raw / sum(fr$p * raw)
  } else rep(1, n_frailty)

  # Per-bin transmissibility frailty, rank-correlated with sus. If sd_trans
  # matches the population sd_pop the bins already sit at the right values;
  # otherwise map ranks through the matching Beta inverse CDF.
  trans_frailty <- if (sd_trans > 0) {
    raw <- if (sd_trans == sd_pop) {
      exp(frailty_amp * fr$x)
    } else {
      cf_pop        <- (0.25 / sd_pop^2)   - 1
      cf_t          <- (0.25 / sd_trans^2) - 1
      ranks         <- pbeta(fr$x, 0.5 * cf_pop, 0.5 * cf_pop)
      exp(frailty_amp * qbeta(ranks, 0.5 * cf_t, 0.5 * cf_t))
    }
    raw / sum(fr$p * raw)
  } else {
    rep(1, n_frailty)
  }

  n_groups       <- 2 * n_frailty
  N_groups       <- c(n_total - vac_counts, vac_counts)
  susceptibility <- c(frailty, alpha * frailty)
  transmissibility <- c(trans_frailty, trans_frailty)  # vaccine acts on sus only
  # Unit mixing entries (NOT /n_groups). The dust SIR model uses
  #   foi = beta * sum_j mm[i,j] * I_j / N
  # so this gives foi = beta * I_total/N in the homogeneous limit, i.e.
  # R0 = beta/gamma — matching the homog SIR wrapper's mm = matrix(1,2,2).
  mm             <- matrix(1, nrow=n_groups, ncol=n_groups)
  
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
                      dt=dt, timepoints=timepoints, n_sim=n_sim, cores=cores,
                      seed=seed)
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

# ---------------------------------------------------------------------------
# Heterogeneous vaccine effect: alpha varies BY INDIVIDUAL
# ---------------------------------------------------------------------------
#
# The frailty models put heterogeneity in baseline susceptibility, expressed
# by everyone. Here baseline susceptibility is 1 for all, and the vaccine's
# effect alpha_i varies between people: it is a fixed individual attribute
# (how well *this* person responds if vaccinated), so it is defined for the
# unvaccinated too and only EXPRESSED by those actually vaccinated. That is
# what makes the flip counterfactual ("what if i had been vaccinated?")
# well-posed -- it needs i's own alpha_i either way.
#
# Population structure is therefore the same 2 x K group layout the frailty
# models use, but with the heterogeneity on the vaccinated side only:
#
#     group k         (k = 1..K)  unvaccinated, bin k, susceptibility 1
#     group K + k                 vaccinated,   bin k, susceptibility alpha_k
#
# A "bin" is a slice of the alpha distribution, a fixed partition of the
# population. Allocation draws which members of each bin get vaccinated
# (multivariate hypergeometric), so varying the allocation varies WHICH
# alpha_i values are expressed -- exactly the sampling variability of
# interest.

# Equal-probability bins of the per-individual vaccine effect.
#
#   alpha_i ~ Beta(mean = alpha, sd = kappa * sqrt(alpha * (1 - alpha)))
#
# Spread is given as a FRACTION kappa of the largest sd a Beta with that mean
# admits, never as an absolute sd. The maximum, sqrt(mu(1-mu)), moves with the
# mean, so an absolute sd turns infeasible the moment the optimiser walks
# alpha towards 0 or 1 -- mid-fit, silently. In this parameterisation the
# concentration is
#       c = mu(1-mu)/s^2 - 1 = 1/kappa^2 - 1,
# free of mu, so any kappa in [0,1) is valid at every alpha. kappa = 0
# reproduces the homogeneous model exactly.
#
# Bins carry equal probability (p_k = 1/K) and sit at the CONDITIONAL MEAN of
# the distribution over each quantile slice, not at the slice midpoint. That
# makes sum_k p_k alpha_k == alpha exactly, by the partial-first-moment
# identity E[X 1{a<X<b}] = mu (I_b(s1+1,s2) - I_a(s1+1,s2)); midpoint bins are
# only mean-accurate to O(K^-2), and that discretisation error would bias the
# fitted alpha rather than just blur it.
#
# Note what does and does not depend on alpha: the bin WEIGHTS are fixed at
# 1/K, only the bin VALUES move. The allocation (who sits in which bin, and
# who is vaccinated) is therefore frozen across the whole fit, and alpha only
# slides the values -- which is what keeps the loss surface smooth enough for
# Nelder-Mead and the finite-difference posterior covariance.
alpha_ve_bins <- function(alpha, kappa = 0, n_alpha = 20L) {
  n_alpha <- as.integer(n_alpha)
  if (n_alpha < 1L) stop("n_alpha must be >= 1")
  if (!is.finite(kappa) || kappa < 0 || kappa >= 1)
    stop("alpha_kappa must be in [0, 1)")
  p <- rep(1 / n_alpha, n_alpha)
  # Degenerate cases: no spread, a single bin, or a mean on the boundary
  # (Beta is a point mass at 0 or 1 there).
  if (kappa == 0 || n_alpha == 1L || alpha <= 0 || alpha >= 1)
    return(list(x = rep(min(max(alpha, 0), 1), n_alpha), p = p))

  cc <- 1 / kappa^2 - 1
  s1 <- alpha * cc
  s2 <- (1 - alpha) * cc
  u  <- seq(0, 1, length.out = n_alpha + 1L)

  # Evaluate each bin edge from whichever tail it is nearer. Once cc < 1 the
  # Beta is U-shaped, and for a small shape parameter the all-lower-tail forms
  # saturate: qbeta() returns exactly 1 for the upper edges and the pbeta()
  # differences collapse to zero, which drags the discretised mean silently
  # away from alpha -- at alpha = 0.95, kappa = 0.95 it lands on 0.25. The
  # failure is driven by min(s1, s2) going small, so it bites hardest exactly
  # where the science is most interesting (a strongly protective vaccine with
  # a wide spread of individual response).
  lo     <- u <= 0.5
  q      <- numeric(n_alpha + 1L)
  q[lo]  <- qbeta(u[lo],      s1, s2)
  q[!lo] <- qbeta(1 - u[!lo], s1, s2, lower.tail = FALSE)
  Fl <- pbeta(q, s1 + 1, s2)                       # lower-tail partial moment
  Fu <- pbeta(q, s1 + 1, s2, lower.tail = FALSE)   # upper-tail partial moment
  d  <- ifelse(u[-1] <= 0.5,
               Fl[-1] - Fl[-(n_alpha + 1L)],
               Fu[-(n_alpha + 1L)] - Fu[-1])
  x  <- alpha * d * n_alpha                        # E[alpha | bin k]
  if (!all(is.finite(x)))
    stop(sprintf("alpha_ve_bins: non-finite bins at alpha = %.6g, kappa = %.4g",
                 alpha, kappa))

  # Centre and shrink, rather than trusting the quantile arithmetic. The mean
  # is the quantity the optimiser identifies, so it has to be exact -- a
  # drifting mean does not look like a bug, it looks like a converged fit at
  # the wrong alpha. Writing the bins as
  #     x_k = alpha + c * (x_k - xbar)
  # makes sum(p * x_k) == alpha by algebra, whatever the tails did, and c is
  # taken as large as keeps every bin inside [0, 1].
  #
  # c < 1 only bites when alpha sits near a boundary, and there it is not a
  # numerical fudge but the real constraint: a variable on [0,1] with mean
  # alpha cannot have sd above sqrt(alpha(1-alpha)), so the requested spread
  # is genuinely unattainable and the honest answer is the widest that fits.
  # Clamping instead would silently move the mean; shrinking preserves it and
  # keeps the shape.
  xbar <- sum(p * x)
  dev  <- x - xbar
  cmax <- 1
  pos  <- dev > 0
  neg  <- dev < 0
  if (any(pos)) cmax <- min(cmax, min((1 - alpha) / dev[pos]))
  if (any(neg)) cmax <- min(cmax, min(alpha / (-dev[neg])))
  cmax <- max(min(cmax, 1), 0)
  x    <- alpha + cmax * dev
  # Guard the floating-point edges only; the mean is exact by construction.
  list(x = pmin(pmax(x, 0), 1), p = p)
}

# Split `total` individuals into bins of probability `p` so the parts sum to
# `total` exactly (largest-remainder). round(total * p) does not: it can be
# off by several people, which quietly changes the population size the fit is
# matching against.
.ve_bin_sizes <- function(total, p) {
  target <- total * p
  n      <- floor(target)
  rem    <- as.integer(round(total - sum(n)))
  if (rem > 0L) {
    for (idx in order(target - n, decreasing = TRUE)[seq_len(rem)])
      n[idx] <- n[idx] + 1L
  }
  as.integer(n)
}

# Distribute I_ini_total index cases over the unvac/vac group vectors,
# proportional to population, keeping the arm split true to the arm sizes.
# Mirrors the seeding in run_stoch_frailty_cd.
#
# NOT .spread_seeds: that name is already taken further down the file by a
# two-argument (total, sizes) helper, and since that definition is sourced
# LAST it silently wins. The collision surfaced as "unused argument".
.spread_arm_seeds <- function(I_ini_total, N_unvac_grp, N_vac_grp) {
  spread <- function(total, sizes) {
    if (total == 0L) return(integer(length(sizes)))
    target <- total * sizes / sum(sizes)
    ini    <- pmin(floor(target), sizes)
    rem    <- total - sum(ini)
    if (rem > 0L) {
      for (idx in order(target - ini, decreasing = TRUE)) {
        if (rem == 0L) break
        if (ini[idx] < sizes[idx]) { ini[idx] <- ini[idx] + 1L; rem <- rem - 1L }
      }
    }
    as.integer(ini)
  }
  total_unvac <- sum(N_unvac_grp); total_vac <- sum(N_vac_grp)
  if (I_ini_total > total_unvac + total_vac)
    stop(sprintf("I_ini_total (%d) exceeds total population (%d).",
                 I_ini_total, total_unvac + total_vac))
  ut          <- I_ini_total * total_unvac / (total_unvac + total_vac)
  unvac_seeds <- min(floor(ut), total_unvac)
  vac_seeds   <- I_ini_total - unvac_seeds
  if (vac_seeds > total_vac) {
    deficit     <- vac_seeds - total_vac
    vac_seeds   <- as.integer(total_vac)
    unvac_seeds <- unvac_seeds + deficit
  } else if (ut - unvac_seeds > 0.5 && unvac_seeds < total_unvac && vac_seeds > 0L) {
    unvac_seeds <- unvac_seeds + 1L; vac_seeds <- vac_seeds - 1L
  }
  c(spread(as.integer(unvac_seeds), N_unvac_grp),
    spread(as.integer(vac_seeds),   N_vac_grp))
}

# Homogeneous-mixing SIR with a per-individual vaccine effect. Sibling of
# run_stoch_frailty_cd; same output columns (time, sim, C1, C2, vac, unvac,
# CRR) so it drops into the same fit machinery.
#
#   N            : per-arm population (total = 2N), matching the frailty
#                  wrapper's convention.
#   alpha_kappa  : spread of alpha_i, as a fraction of the maximum (see
#                  alpha_ve_bins). 0 = homogeneous, identical to the plain
#                  SIR wrapper.
#   vac_counts   : vaccinated count per alpha bin. Pass the allocation drawn
#                  once per config; NULL gives the deterministic round(f * n_k)
#                  split, i.e. no allocation variability.
run_stoch_ve_hetero_cd <- function(beta = 1, f = 0.5, N = 1000, t = 100,
                                   alpha_kappa = 0, n_alpha = 20L,
                                   gamma = 1 / 2, vac_counts = NULL,
                                   I_ini_total = 1, timepoints = seq(0, t, 1),
                                   n_sim = 100, cores = 10, dt = 0.1,
                                   susceptibility = c(1, 0.5), seed = NULL) {
  alpha <- susceptibility[2]
  ab    <- alpha_ve_bins(alpha, alpha_kappa, n_alpha)
  K     <- length(ab$x)

  n_total_k <- .ve_bin_sizes(round(2 * N), ab$p)
  if (is.null(vac_counts)) vac_counts <- round(f * n_total_k)
  vac_counts <- pmin(pmax(as.integer(vac_counts), 0L), n_total_k)

  n_groups <- 2L * K
  N_groups <- c(n_total_k - vac_counts, vac_counts)
  # Baseline susceptibility for everyone; the vaccinated express alpha_k.
  sus      <- c(rep(susceptibility[1], K), ab$x)
  trans    <- rep(1, n_groups)         # vaccine acts on susceptibility only
  # Unit mixing entries (NOT /n_groups): the dust SIR model uses
  # foi = beta * sum_j mm[i,j] * I_j / N, so this keeps R0 = beta/gamma in the
  # homogeneous limit -- same convention as run_stoch_frailty_cd.
  mm       <- matrix(1, nrow = n_groups, ncol = n_groups)

  I_ini <- .spread_arm_seeds(I_ini_total,
                         N_groups[seq_len(K)],
                         N_groups[(K + 1L):n_groups])

  raw <- run_stoch_cd_dust(mm, beta = beta, N = N_groups, t = t, I_ini = I_ini,
                           susceptibility = sus, transmissibility = trans,
                           gamma = gamma, dt = dt, timepoints = timepoints,
                           n_sim = n_sim, cores = cores, seed = seed)
  setDT(raw)
  N_vac_tot   <- sum(vac_counts)
  N_unvac_tot <- sum(n_total_k - vac_counts)
  raw[, vac   := rowSums(.SD), .SDcols = paste0("C", (K + 1L):n_groups)]
  raw[, unvac := rowSums(.SD), .SDcols = paste0("C", 1:K)]
  raw[, CRR   := (vac / N_vac_tot) / (unvac / N_unvac_tot)]
  raw[, C1 := unvac]
  raw[, C2 := vac]
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


# ---------------------------------------------------------------------------
# Stochastic EATE via frozen-field counterfactual (network)
# ---------------------------------------------------------------------------
#
# Replaces the N+1 deterministic perturbation runs of get_eate_network
# with a stochastic factual averaged over n_rep replicates plus a
# per-individual frozen-field counterfactual. Each EATE expectation
# still uses 2N inputs but only N of them come from new "simulations":
# the N counterfactual expectations (the ones that flip a unit's vac
# status) are derived analytically from extracted force-of-infection
# trajectories.
#
#   For each replicate r and each individual i accumulate the FOI they
#   would have experienced assuming they stayed susceptible:
#       cumFOI_i^(r)(t) = integral_0^t  beta/k_mean * sum_{j in N(i)} I_j^(r)(s) ds
#   Then under counterfactual susceptibility sus_v,
#       P_i^v(t) = 1 - mean_r exp( - sus_v * cumFOI_i^(r)(t) )
#   with sus_v = 1 ("unvac counterfactual") or alpha ("vac counterfactual",
#   matching the codebase convention that susceptibility[2] = alpha).
#
#   Like the deterministic frozen-field EATE, individuals whose factual
#   vac status matches the counterfactual contribute their *factual*
#   case probability (averaged over reps), and the flipped side uses the
#   frozen counterfactual:
#       num   = sum_{i in vac}     P_i^factual  +  sum_{i not in vac} P_i^vac
#       denom = sum_{i not in vac} P_i^factual  +  sum_{i in vac}     P_i^unvac
#       EATE  = num / denom
#   CRR(t) is taken from the factual cumulative cases per group.
#
# Outer parallelism: parallel::mclapply across n_vac allocations
# (mc.cores). Inner parallelism: dust2 threads inside one replicate
# batch (inner_cores). Defaults send all cores to the outer loop.

# Helper: reshape one column of a (time outer, sim inner) data.table
# into an [n_t, n_rep] matrix. Mirrors the dust2 output layout used by
# run_stoch_adj / run_stoch_cd_dust.
.dt_col_to_t_rep_matrix <- function(v, n_t, n_rep) {
  matrix(v, nrow = n_t, ncol = n_rep, byrow = TRUE)
}

# Trapezoidal cumulative integral along the first axis of a [n_t, n_x]
# matrix `FI` against time grid `timepoints`. Returns [n_t, n_x] with
# row 1 = 0.
.cum_trapz <- function(FI, timepoints) {
  n_t <- nrow(FI)
  out <- matrix(0, nrow = n_t, ncol = ncol(FI))
  if (n_t < 2L) return(out)
  dt_vec <- diff(timepoints)
  for (it in seq.int(2L, n_t)) {
    out[it, ] <- out[it - 1L, ] + (FI[it, ] + FI[it - 1L, ]) / 2 * dt_vec[it - 1L]
  }
  out
}

get_stoch_eate_network <- function(beta = 1, susceptibility = c(1, 1), f = 0.5,
                                   N = 200, t = 15, pl_alpha = 3, c_ij = NULL,
                                   n_vac = 10, n_rep = 20,
                                   k_mean = 6, gamma = 1 / 3, dt = 0.1,
                                   timepoints = NULL, init_I = 2,
                                   mc.cores = 10, inner_cores = 1,
                                   vac_list = NULL, adj = NULL) {
  alpha <- susceptibility[2]
  if (is.null(c_ij) && is.null(adj))
    c_ij <- get_conact_matrix_pl(N, pl_alpha, mean_k = k_mean)
  # Build the adjacency ONCE, not once per allocation (see run_stoch_adj).
  if (is.null(adj)) adj <- contact_matrix_to_adj(c_ij)
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t <- length(timepoints)

  # If explicit allocations are supplied, use them in order and override
  # n_vac to match; otherwise sample a fresh vac set per iteration.
  if (!is.null(vac_list)) n_vac <- length(vac_list)

  run_one_allocation <- function(i) {
    vac      <- if (!is.null(vac_list)) vac_list[[i]]
                else sample(seq_len(N), round(f * N))
    non_vac  <- setdiff(seq_len(N), vac)
    sim_id   <- runif(1)

    susept       <- rep(1, N)
    susept[vac]  <- alpha
    I_ini_vec    <- c(rep(1L, init_I), rep(0L, N - init_I))

    raw <- run_stoch_adj(c_ij,
                         beta           = N * beta / k_mean,
                         t              = t,
                         I_ini          = I_ini_vec,
                         susceptibility = susept,
                         gamma          = gamma,
                         dt             = dt,
                         timepoints     = timepoints,
                         n_sim          = n_rep,
                         cores          = inner_cores,
                         adj            = adj)
    setDT(raw)

    # Factual cumulative cases per (timepoint, individual), averaged over
    # replicates. C_k is 0/1 per replicate; rowMeans gives P_factual_i.
    P_factual <- matrix(0, nrow = n_t, ncol = N)
    I_mat     <- array(0, dim = c(n_t, n_rep, N))
    for (k in seq_len(N)) {
      Ck <- .dt_col_to_t_rep_matrix(raw[[paste0("C", k)]], n_t, n_rep)
      P_factual[, k] <- rowMeans(Ck)
      I_mat[,, k]    <- .dt_col_to_t_rep_matrix(raw[[paste0("I", k)]], n_t, n_rep)
    }

    # Per-replicate cumulative FOI per individual:
    #   FOI_i^(r)(t) = beta/k_mean * sum_{j in N(i)} I_j^(r)(t)
    cum_foi_traj <- array(0, dim = c(n_t, n_rep, N))
    for (r in seq_len(n_rep)) {
      I_traj_r <- I_mat[, r, ]                                 # [n_t, N]
      FI_r     <- (I_traj_r %*% t(c_ij)) * (beta / k_mean)     # [n_t, N]
      cum_foi_traj[, r, ] <- .cum_trapz(FI_r, timepoints)
    }

    # Hybrid EATE: matching side uses P_factual, flipped side uses the
    # frozen counterfactual averaged over replicates. num / denom are
    # the counterfactual "sum Y_i(vac)" and "sum Y_i(unvac)" totals; VE
    # = 1 - num/denom on the ratio scale, AVE = (denom - num) / N on the
    # per-person absolute scale.
    eate_t <- numeric(n_t); num_t <- numeric(n_t); denom_t <- numeric(n_t)
    ave_t  <- numeric(n_t); crr_t <- numeric(n_t); crr_ave_t <- numeric(n_t)
    for (it in seq_len(n_t)) {
      cfi        <- cum_foi_traj[it, , ]                       # [n_rep, N]
      P_vac_cf   <- 1 - colMeans(exp(-alpha * cfi))            # length N
      P_unvac_cf <- 1 - colMeans(exp(-cfi))                    # length N
      P_fac      <- P_factual[it, ]                            # length N
      num   <- sum(P_fac[vac])     + sum(P_vac_cf[non_vac])
      denom <- sum(P_fac[non_vac]) + sum(P_unvac_cf[vac])
      num_t[it]   <- num
      denom_t[it] <- denom
      eate_t[it]  <- num / denom
      ave_t[it]   <- (denom - num) / N
      ar_fac_vac   <- sum(P_fac[vac])     / length(vac)
      ar_fac_unvac <- sum(P_fac[non_vac]) / length(non_vac)
      crr_t[it]     <- ar_fac_vac / ar_fac_unvac
      crr_ave_t[it] <- ar_fac_unvac - ar_fac_vac
    }

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = crr_t, ave = crr_ave_t,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  res <- parallel::mclapply(seq_len(n_vac),
                            function(i) run_one_allocation(i),
                            mc.cores = mc.cores)
  rbindlist(res, fill = TRUE)
}


# ---------------------------------------------------------------------------
# Stochastic EATE via frozen-field counterfactual (frailty / homogeneous)
# ---------------------------------------------------------------------------
#
# Same hybrid form as get_stoch_eate_network, but the population is split
# into 2*n_frailty groups (n_frailty unvac bins + n_frailty vac bins).
# Mixing is uniform, so the FOI experienced by any individual is identical
# at each (t, replicate) — we accumulate a single cum_foi_rep[t, r]
# trajectory per replicate and broadcast it across bins.
#
#   FOI(t, r) = beta * sum_g trans_g * I_g^(r)(t) / N_total
#   cum_foi_rep[t, r] = trapezoidal integral of FOI over timepoints
#   sus_v_k = frailty[k] for v=0 (unvac), alpha * frailty[k] for v=1 (vac)
#   P_v_k(t) = 1 - mean_r exp( - sus_v_k * cum_foi_rep[t, r] )
#   num   = sum_k [ C_vac_bin[t, k]   + P_vac_k(t)   * N_unvac_grp[k] ]
#   denom = sum_k [ C_unvac_bin[t, k] + P_unvac_k(t) * N_vac_grp[k]   ]
#   EATE  = num / denom
#
# The "no heterogeneity" case is n_frailty = 1 with sd = sd_trans = 0:
# the population has 2 groups (unvac + vac), the formula simplifies to
# the user's e^{-(1 - alpha * v) * s * FOI} form with the codebase
# convention sus_v = (1-v) + alpha * v and s = frailty[k] = 1.

get_stoch_eate_frailty <- function(alpha, sd = 0, sd_trans = 0, beta = 1, R = NULL,
                                   f = 0.5, N = 1000, t = 30, n_frailty = 1,
                                   n_vac = 10, n_rep = 20,
                                   gamma = 1, dt = 0.1, timepoints = NULL,
                                   I_ini_total = 1, mc.cores = 10,
                                   inner_cores = 1, frailty_amp = 2.5) {
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t <- length(timepoints)

  if (n_frailty == 1L) {
    if (sd > 0 || sd_trans > 0) {
      warning("n_frailty=1 collapses bins; sd / sd_trans inputs are ignored.")
    }
    fr_p          <- 1
    frailty       <- 1
    trans_frailty <- 1
  } else {
    sd_pop <- if (sd > 0) sd else sd_trans
    if (sd_pop <= 0) {
      stop("get_stoch_eate_frailty with n_frailty>1 requires sd>0 or sd_trans>0")
    }
    fr      <- get_frailty(sd = sd_pop, n = n_frailty)
    fr_p    <- fr$p
    # Frailty arrays normalised to population-mean 1; see run_stoch_frailty_cd.
    frailty <- if (sd > 0) {
      raw <- exp(frailty_amp * fr$x); raw / sum(fr$p * raw)
    } else rep(1, n_frailty)
    trans_frailty <- if (sd_trans > 0) {
      raw <- if (sd_trans == sd_pop) {
        exp(frailty_amp * fr$x)
      } else {
        cf_pop <- (0.25 / sd_pop^2)   - 1
        cf_t   <- (0.25 / sd_trans^2) - 1
        ranks  <- pbeta(fr$x, 0.5 * cf_pop, 0.5 * cf_pop)
        exp(frailty_amp * qbeta(ranks, 0.5 * cf_t, 0.5 * cf_t))
      }
      raw / sum(fr$p * raw)
    } else {
      rep(1, n_frailty)
    }
  }

  n_total_k     <- round(2 * N * fr_p)
  n_groups      <- 2L * n_frailty
  trans_all     <- c(trans_frailty, trans_frailty)
  sus_unvac_bin <- frailty
  sus_vac_bin   <- alpha * frailty
  N_total       <- sum(n_total_k)

  if (!is.null(R)) {
    beta <- get_beta(R, alpha, sd, sd_trans = sd_trans, f = f,
                     N = N, n_frailty = n_frailty, gamma = gamma)
  }

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

  run_one_allocation <- function() {
    n_vac_total <- round(f * N_total)
    vac_counts  <- if (n_frailty == 1L) {
      as.integer(n_vac_total)
    } else {
      tabulate(sample(rep(seq_len(n_frailty), n_total_k), n_vac_total),
               nbins = n_frailty)
    }
    sim_id <- runif(1)

    N_groups <- c(n_total_k - vac_counts, vac_counts)
    susceptibility   <- c(frailty,       alpha * frailty)
    transmissibility <- c(trans_frailty, trans_frailty)
    # Unit mixing entries (NOT /n_groups) so beta=R0*gamma in the
    # homogeneous limit — matches run_stoch_frailty_cd.
    mm <- matrix(1, nrow = n_groups, ncol = n_groups)

    N_unvac_grp <- N_groups[seq_len(n_frailty)]
    N_vac_grp   <- N_groups[(n_frailty + 1L):n_groups]
    total_unvac <- sum(N_unvac_grp); total_vac <- sum(N_vac_grp)

    ut          <- I_ini_total * total_unvac / (total_unvac + total_vac)
    unvac_seeds <- min(floor(ut), total_unvac)
    vac_seeds   <- I_ini_total - unvac_seeds
    if (vac_seeds > total_vac) {
      deficit     <- vac_seeds - total_vac
      vac_seeds   <- as.integer(total_vac)
      unvac_seeds <- unvac_seeds + deficit
    } else if (ut - unvac_seeds > 0.5 && unvac_seeds < total_unvac && vac_seeds > 0L) {
      unvac_seeds <- unvac_seeds + 1L; vac_seeds <- vac_seeds - 1L
    }
    I_ini <- c(spread(as.integer(unvac_seeds), N_unvac_grp),
               spread(as.integer(vac_seeds),   N_vac_grp))

    raw <- run_stoch_cd_dust(mm, beta = beta, N = N_groups, t = t,
                             I_ini = I_ini,
                             susceptibility = susceptibility,
                             transmissibility = transmissibility,
                             gamma = gamma, dt = dt,
                             timepoints = timepoints,
                             n_sim = n_rep, cores = inner_cores)
    setDT(raw)

    # I_mat[t, r, g] for all 2*n_frailty groups.
    I_mat <- array(0, dim = c(n_t, n_rep, n_groups))
    for (g in seq_len(n_groups)) {
      I_mat[,, g] <- .dt_col_to_t_rep_matrix(raw[[paste0("I", g)]], n_t, n_rep)
    }

    # Uniform mixing => one common FOI per (t, rep). Matches the
    # deterministic frozen-field normalisation in get_frailty_eate.
    cum_foi_rep <- matrix(0, nrow = n_t, ncol = n_rep)
    for (r in seq_len(n_rep)) {
      I_traj_r     <- if (n_groups == 1L) {
        matrix(I_mat[, r, ], nrow = n_t, ncol = 1)
      } else {
        I_mat[, r, ]                                              # [n_t, n_groups]
      }
      I_weighted_r <- as.numeric(I_traj_r %*% trans_all)          # [n_t]
      FI_r         <- beta * I_weighted_r / N_total               # [n_t]
      cum_foi_rep[, r] <- .cum_trapz(matrix(FI_r, ncol = 1), timepoints)[, 1]
    }

    # Factual cumulative cases per bin, averaged over reps.
    C_unvac_bin <- matrix(0, nrow = n_t, ncol = n_frailty)
    C_vac_bin   <- matrix(0, nrow = n_t, ncol = n_frailty)
    for (k in seq_len(n_frailty)) {
      C_unvac_bin[, k] <- rowMeans(.dt_col_to_t_rep_matrix(
        raw[[paste0("C", k)]], n_t, n_rep))
      C_vac_bin[, k]   <- rowMeans(.dt_col_to_t_rep_matrix(
        raw[[paste0("C", n_frailty + k)]], n_t, n_rep))
    }

    # Hybrid EATE per timepoint: matching side uses factual cases, flipped
    # side uses the per-bin frozen counterfactual scaled by bin population.
    # Also compute AVE = (denom - num) / N_total (per-person absolute effect).
    eate_t <- numeric(n_t); num_t <- numeric(n_t); denom_t <- numeric(n_t)
    ave_t  <- numeric(n_t)
    for (it in seq_len(n_t)) {
      cfi <- cum_foi_rep[it, ]                                    # [n_rep]
      # outer(cfi, sus_bin) -> [n_rep, n_frailty]; colMeans gives [n_frailty].
      surv_vac_k   <- colMeans(exp(-outer(cfi, sus_vac_bin)))
      surv_unvac_k <- colMeans(exp(-outer(cfi, sus_unvac_bin)))
      P_vac_k      <- 1 - surv_vac_k
      P_unvac_k    <- 1 - surv_unvac_k

      num   <- sum(C_vac_bin[it, ])   + sum(P_vac_k   * N_unvac_grp)
      denom <- sum(C_unvac_bin[it, ]) + sum(P_unvac_k * N_vac_grp)
      num_t[it]   <- num
      denom_t[it] <- denom
      eate_t[it]  <- num / denom
      ave_t[it]   <- (denom - num) / N_total
    }
    ar_fac_vac_t   <- rowSums(C_vac_bin)   / total_vac
    ar_fac_unvac_t <- rowSums(C_unvac_bin) / total_unvac
    crr_t     <- ar_fac_vac_t / ar_fac_unvac_t
    crr_ave_t <- ar_fac_unvac_t - ar_fac_vac_t

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = crr_t, ave = crr_ave_t,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  res <- parallel::mclapply(seq_len(n_vac),
                            function(i) run_one_allocation(),
                            mc.cores = mc.cores)
  rbindlist(res, fill = TRUE)
}


# ---------------------------------------------------------------------------
# Stochastic EATE with a per-individual vaccine effect
# ---------------------------------------------------------------------------
#
# Frozen-field counterfactual for run_stoch_ve_hetero_cd, built like
# get_stoch_eate_frailty: run one allocation, freeze the realised cumulative
# force of infection, then ask what would have happened to each bin under the
# flipped treatment.
#
# The one thing that differs from every other EATE here: the counterfactual
# is per BIN, because alpha_i is an individual attribute. Flipping an
# unvaccinated person in bin k to vaccinated gives them THEIR OWN alpha_k, not
# a population-average alpha:
#
#     P(case | vaccinated,   bin k) = 1 - E_reps[ exp(-alpha_k * cumFOI) ]
#     P(case | unvaccinated, bin k) = 1 - E_reps[ exp(-1       * cumFOI) ]
#
# so the vaccinated-side counterfactual varies across bins even though the
# frozen field does not. Averaging alpha first and exponentiating after would
# be wrong by Jensen's inequality -- exp(-a x) is convex in a, so a single
# mean alpha understates the counterfactual case count, and the size of that
# error grows with alpha_kappa. That is precisely the effect being measured,
# so getting it from the average would hide the whole phenomenon.
#
# Allocation variability is the outer loop (n_vac draws): each draw
# re-samples which members of each bin are vaccinated, so it samples which
# alpha_i values are expressed. The alpha_i values themselves are a fixed
# population attribute and do NOT get redrawn.
get_stoch_eate_ve_hetero <- function(alpha, alpha_kappa = 0, n_alpha = 20L,
                                     beta = 1, f = 0.5, N = 1000, t = 30,
                                     n_vac = 10, n_rep = 20,
                                     gamma = 1, dt = 0.1, timepoints = NULL,
                                     I_ini_total = 1, mc.cores = 10,
                                     inner_cores = 1) {
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t <- length(timepoints)

  ab  <- alpha_ve_bins(alpha, alpha_kappa, n_alpha)
  K   <- length(ab$x)
  n_total_k <- .ve_bin_sizes(round(2 * N), ab$p)
  N_total   <- sum(n_total_k)
  n_groups  <- 2L * K

  sus_unvac_bin <- rep(1, K)      # baseline susceptibility, everyone
  sus_vac_bin   <- ab$x           # this bin's own vaccine effect
  trans_all     <- rep(1, n_groups)

  run_one_allocation <- function() {
    n_vac_total <- round(f * N_total)
    # Multivariate hypergeometric: which members of each bin are vaccinated.
    vac_counts  <- if (K == 1L) {
      as.integer(n_vac_total)
    } else {
      tabulate(sample(rep(seq_len(K), n_total_k), n_vac_total), nbins = K)
    }
    sim_id <- runif(1)

    N_groups    <- c(n_total_k - vac_counts, vac_counts)
    N_unvac_grp <- N_groups[seq_len(K)]
    N_vac_grp   <- N_groups[(K + 1L):n_groups]
    total_unvac <- sum(N_unvac_grp); total_vac <- sum(N_vac_grp)

    mm    <- matrix(1, nrow = n_groups, ncol = n_groups)
    I_ini <- .spread_arm_seeds(I_ini_total, N_unvac_grp, N_vac_grp)

    raw <- run_stoch_cd_dust(mm, beta = beta, N = N_groups, t = t,
                             I_ini = I_ini,
                             susceptibility = c(sus_unvac_bin, sus_vac_bin),
                             transmissibility = trans_all,
                             gamma = gamma, dt = dt, timepoints = timepoints,
                             n_sim = n_rep, cores = inner_cores)
    setDT(raw)

    I_mat <- array(0, dim = c(n_t, n_rep, n_groups))
    for (g in seq_len(n_groups))
      I_mat[,, g] <- .dt_col_to_t_rep_matrix(raw[[paste0("I", g)]], n_t, n_rep)

    # Uniform mixing => one common FOI per (t, rep).
    cum_foi_rep <- matrix(0, nrow = n_t, ncol = n_rep)
    for (r in seq_len(n_rep)) {
      I_traj_r <- if (n_groups == 1L) matrix(I_mat[, r, ], nrow = n_t, ncol = 1)
                  else                I_mat[, r, ]
      FI_r     <- beta * as.numeric(I_traj_r %*% trans_all) / N_total
      cum_foi_rep[, r] <- .cum_trapz(matrix(FI_r, ncol = 1), timepoints)[, 1]
    }

    C_unvac_bin <- matrix(0, nrow = n_t, ncol = K)
    C_vac_bin   <- matrix(0, nrow = n_t, ncol = K)
    for (k in seq_len(K)) {
      C_unvac_bin[, k] <- rowMeans(.dt_col_to_t_rep_matrix(
        raw[[paste0("C", k)]], n_t, n_rep))
      C_vac_bin[, k]   <- rowMeans(.dt_col_to_t_rep_matrix(
        raw[[paste0("C", K + k)]], n_t, n_rep))
    }

    eate_t <- numeric(n_t); num_t <- numeric(n_t); denom_t <- numeric(n_t)
    ave_t  <- numeric(n_t)
    for (it in seq_len(n_t)) {
      cfi <- cum_foi_rep[it, ]                               # [n_rep]
      # outer(cfi, sus_bin) -> [n_rep, K]; colMeans averages over reps only,
      # keeping the per-bin alpha_k intact.
      P_vac_k   <- 1 - colMeans(exp(-outer(cfi, sus_vac_bin)))
      P_unvac_k <- 1 - colMeans(exp(-outer(cfi, sus_unvac_bin)))

      num   <- sum(C_vac_bin[it, ])   + sum(P_vac_k   * N_unvac_grp)
      denom <- sum(C_unvac_bin[it, ]) + sum(P_unvac_k * N_vac_grp)
      num_t[it]   <- num
      denom_t[it] <- denom
      eate_t[it]  <- num / denom
      ave_t[it]   <- (denom - num) / N_total
    }
    ar_fac_vac_t   <- rowSums(C_vac_bin)   / total_vac
    ar_fac_unvac_t <- rowSums(C_unvac_bin) / total_unvac
    crr_t     <- ar_fac_vac_t / ar_fac_unvac_t
    crr_ave_t <- ar_fac_unvac_t - ar_fac_vac_t

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = crr_t, ave = crr_ave_t,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  res <- parallel::mclapply(seq_len(n_vac),
                            function(i) run_one_allocation(),
                            mc.cores = mc.cores)
  rbindlist(res, fill = TRUE)
}


# ---------------------------------------------------------------------------
# Stochastic EATE via frozen-field counterfactual (linear)
# ---------------------------------------------------------------------------
#
# Linear (exposure-only) model: each individual is infected independently
# at rate beta * susceptibility, no FOI from infected. The cumulative FOI
# experienced is therefore deterministic — cum_foi(t) = beta * t — and the
# counterfactual case probabilities are analytical:
#   P_v(t) = 1 - exp(-sus_v * beta * t)
# Stochastic noise enters only via the factual case counts per replicate.
# The hybrid EATE mirrors get_stoch_eate_network: factually matching
# individuals contribute their factual case probability (mean over reps),
# flipped individuals contribute the frozen counterfactual.
get_stoch_eate_linear <- function(beta = 1, susceptibility = c(1, 1), f = 0.5,
                                  N = 200, t = 30,
                                  n_vac = 10, n_rep = 20,
                                  dt = 0.1, timepoints = NULL,
                                  mc.cores = 10, inner_cores = 1, seed = NULL) {
  alpha <- susceptibility[2]
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t     <- length(timepoints)
  N_unvac <- round(N * (1 - f))
  N_vac   <- N - N_unvac

  # cumFOI is deterministic for the linear model.
  cum_foi    <- beta * timepoints
  P_unvac_cf <- 1 - exp(-1     * cum_foi)
  P_vac_cf   <- 1 - exp(-alpha * cum_foi)

  run_one_allocation <- function() {
    sim_id <- runif(1)
    raw <- run_stoch_linear_dust(
      beta = beta, N = c(N_unvac, N_vac),
      susceptibility = c(1, alpha),
      t = t, dt = dt, timepoints = timepoints,
      n_sim = n_rep, cores = inner_cores, seed = seed)
    setDT(raw)

    C1_mat <- .dt_col_to_t_rep_matrix(raw$C1, n_t, n_rep)
    C2_mat <- .dt_col_to_t_rep_matrix(raw$C2, n_t, n_rep)
    P_fac_unvac <- rowMeans(C1_mat) / N_unvac
    P_fac_vac   <- rowMeans(C2_mat) / N_vac

    num_t   <- N_vac   * P_fac_vac   + N_unvac * P_vac_cf
    denom_t <- N_unvac * P_fac_unvac + N_vac   * P_unvac_cf
    eate_t  <- num_t / denom_t
    ave_t   <- (denom_t - num_t) / N
    crr_t     <- P_fac_vac / P_fac_unvac
    crr_ave_t <- P_fac_unvac - P_fac_vac

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = crr_t, ave = crr_ave_t,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  res <- parallel::mclapply(seq_len(n_vac),
                            function(i) run_one_allocation(),
                            mc.cores = mc.cores)
  rbindlist(res, fill = TRUE)
}


# ---------------------------------------------------------------------------
# Stochastic EATE via frozen-field counterfactual (homogeneous SIR)
# ---------------------------------------------------------------------------
#
# Two-group SIR (vac + unvac, mass-action mixing). FOI on every
# susceptible is identical at each (t, replicate):
#   FOI^(r)(t) = beta * (I_unvac^(r)(t) + I_vac^(r)(t)) / N
# Cumulative FOI is integrated by trapezoidal rule across timepoints, and
# the hybrid EATE mirrors get_stoch_eate_network / get_stoch_eate_frailty:
#   num   = N_vac   * P_fac_vac   + N_unvac * P_vac_cf
#   denom = N_unvac * P_fac_unvac + N_vac   * P_unvac_cf
#   EATE  = num / denom
# Matches get_EATEs frozen-field formula but with the deterministic
# factual full_0 replaced by an n_rep-average factual stochastic run.
get_stoch_eate_sir <- function(beta = 1, susceptibility = c(1, 1), f = 0.5,
                               N = 200, t = 30, gamma = 1, I_ini = c(2, 2),
                               n_vac = 10, n_rep = 20,
                               dt = 0.1, timepoints = NULL,
                               mc.cores = 10, inner_cores = 1, seed = NULL) {
  alpha <- susceptibility[2]
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t     <- length(timepoints)
  N_unvac <- round(N * (1 - f))
  N_vac   <- N - N_unvac

  run_one_allocation <- function() {
    sim_id <- runif(1)
    raw <- run_stoch_cd_dust(
      matrix(rep(1, 4), nrow = 2),
      beta = beta, N = c(N_unvac, N_vac),
      t = t, I_ini = I_ini,
      susceptibility = c(1, alpha),
      gamma = gamma, dt = dt,
      timepoints = timepoints,
      n_sim = n_rep, cores = inner_cores, seed = seed)
    setDT(raw)

    I1_mat      <- .dt_col_to_t_rep_matrix(raw$I1, n_t, n_rep)
    I2_mat      <- .dt_col_to_t_rep_matrix(raw$I2, n_t, n_rep)
    I_total_mat <- I1_mat + I2_mat
    FI_mat      <- beta * I_total_mat / N
    cum_foi     <- .cum_trapz(FI_mat, timepoints)

    P_vac_cf   <- 1 - rowMeans(exp(-alpha * cum_foi))
    P_unvac_cf <- 1 - rowMeans(exp(-1     * cum_foi))

    C1_mat <- .dt_col_to_t_rep_matrix(raw$C1, n_t, n_rep)
    C2_mat <- .dt_col_to_t_rep_matrix(raw$C2, n_t, n_rep)
    P_fac_unvac <- rowMeans(C1_mat) / N_unvac
    P_fac_vac   <- rowMeans(C2_mat) / N_vac

    num_t   <- N_vac   * P_fac_vac   + N_unvac * P_vac_cf
    denom_t <- N_unvac * P_fac_unvac + N_vac   * P_unvac_cf
    eate_t  <- num_t / denom_t
    ave_t   <- (denom_t - num_t) / N
    crr_t     <- P_fac_vac / P_fac_unvac
    crr_ave_t <- P_fac_unvac - P_fac_vac

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = crr_t, ave = crr_ave_t,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  res <- parallel::mclapply(seq_len(n_vac),
                            function(i) run_one_allocation(),
                            mc.cores = mc.cores)
  rbindlist(res, fill = TRUE)
}

# ---------------------------------------------------------------------------
# Multi-site RCT helpers (shared by the fit simulator in run_fit_array and the
# EATE below). L sites, each a full randomised 2-group trial population laid
# out block-diagonally so sites share no force-of-infection. Groups are
# interleaved [unvac_1, vac_1, unvac_2, vac_2, ...]; the within-site mixing
# weight = sum(N)/N_site = L so each site keeps R0 = beta/gamma.
#
# site_icc in [0,1] is the intra-site correlation of vaccine status. Per-site
# counts are Beta-binomial: a per-site probability p_l ~ Beta(f*kappa,
# (1-f)*kappa) with kappa = (1-icc)/icc, then Binomial(N_site, p_l), nudged to
# the exact global total N_vac.
#   icc = 0 -> p_l = f, counts ~ Binomial(N_site, f): SIMPLE (individual)
#             randomisation. Per-site counts fluctuate by chance around f with
#             the global total held exact; the between-allocation spread at
#             icc=0 is this chance imbalance. Design effect ~ 1.
#   0<icc<1 -> extra site-level over-dispersion (partial clustering).
#   icc = 1 -> whole sites are single-arm (cluster randomised).
# ---------------------------------------------------------------------------

# Integer partition of N_tot into L (near-equal) site sizes.
multisite_site_sizes <- function(N_tot, L) {
  as.integer(diff(round(seq(0, N_tot, length.out = L + 1L))))
}

# Per-site vaccinated counts. Uses the current RNG stream, so set the seed in
# the caller for reproducibility. Returns integers in [0, N_site_vec] summing
# exactly to N_vac_total.
multisite_vac_counts <- function(N_site_vec, N_vac_total, f, icc) {
  L <- length(N_site_vec)
  p <- if (icc <= 0) {
    rep(f, L)
  } else if (icc >= 1) {
    n_vac_sites <- round(f * L)
    pv <- rep(0, L)
    if (n_vac_sites > 0L) pv[sample.int(L, n_vac_sites)] <- 1
    pv
  } else {
    kappa <- (1 - icc) / icc
    rbeta(L, f * kappa, (1 - f) * kappa)
  }
  # Random per-site counts (Binomial), so icc=0 is simple randomisation with
  # per-site chance imbalance rather than deterministic exact balance.
  vac  <- rbinom(L, N_site_vec, pmin(pmax(p, 0), 1))
  # Nudge one unit at a time to hit N_vac_total exactly, respecting caps.
  step  <- sign(N_vac_total - sum(vac))
  guard <- 0L
  while (sum(vac) != N_vac_total && guard < 1000L * L) {
    elig <- if (step > 0) which(vac < N_site_vec) else which(vac > 0L)
    if (!length(elig)) break
    i <- elig[sample.int(length(elig), 1L)]
    vac[i] <- vac[i] + step
    guard  <- guard + 1L
  }
  as.integer(vac)
}

# Distribute `total` initial infecteds across groups proportional to size.
.spread_seeds <- function(total, sizes) {
  if (total <= 0L || sum(sizes) == 0L) return(integer(length(sizes)))
  target <- total * sizes / sum(sizes)
  ini    <- pmin(floor(target), sizes)
  rem    <- total - sum(ini)
  if (rem > 0L) for (idx in order(target - ini, decreasing = TRUE)) {
    if (rem == 0L) break
    if (ini[idx] < sizes[idx]) { ini[idx] <- ini[idx] + 1L; rem <- rem - 1L }
  }
  as.integer(ini)
}

# 2L x 2L block-diagonal contact matrix, interleaved unvac/vac per site.
multisite_block_matrix <- function(L, w) {
  n  <- 2L * L
  mm <- matrix(0, n, n)
  for (l in seq_len(L)) {
    idx <- c(2L * l - 1L, 2L * l)
    mm[idx, idx] <- w
  }
  mm
}

# Multi-site RCT EATE. L independent sites (block-diagonal), each randomised
# with per-site vaccine fraction set by site_icc; factual cases pooled across
# sites, counterfactuals use each site's own FOI. icc=0 is the individually
# randomised design (within-site cancellation); icc=1 is the cluster /
# fully-separated extreme. Structure mirrors get_stoch_eate_sir per site.
get_stoch_eate_sir_multisite <- function(beta = 1, susceptibility = c(1, 1),
                                         f = 0.5, N = 200, t = 30, gamma = 1,
                                         I_ini = c(2, 2), n_sites = 2,
                                         site_icc = 0, n_vac = 10, n_rep = 20,
                                         dt = 0.1, timepoints = NULL,
                                         mc.cores = 10, inner_cores = 1,
                                         seed = NULL) {
  alpha <- susceptibility[2]
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t         <- length(timepoints)
  N_tot       <- round(N)
  L           <- n_sites
  N_vac_total <- round(N_tot * f)
  N_site_vec  <- multisite_site_sizes(N_tot, L)
  mm          <- multisite_block_matrix(L, L)   # within-site weight = L
  I_total     <- sum(I_ini)

  run_one_allocation <- function() {
    sim_id  <- runif(1)
    vac_l   <- multisite_vac_counts(N_site_vec, N_vac_total, f, site_icc)
    unvac_l <- N_site_vec - vac_l
    N_grp   <- as.integer(rbind(unvac_l, vac_l))          # interleaved
    sus_grp <- as.numeric(rbind(rep(1, L), rep(alpha, L)))
    I_grp   <- .spread_seeds(I_total, N_grp)

    raw <- run_stoch_cd_dust(mm, beta = beta, N = N_grp, t = t, I_ini = I_grp,
                             susceptibility = sus_grp, gamma = gamma, dt = dt,
                             timepoints = timepoints, n_sim = n_rep,
                             cores = inner_cores, seed = seed)
    setDT(raw)

    num_t <- numeric(n_t); denom_t <- numeric(n_t)
    tot_Cvac <- numeric(n_t); tot_Cunvac <- numeric(n_t)
    for (l in seq_len(L)) {
      gu <- 2L * l - 1L; gv <- 2L * l
      Iu <- .dt_col_to_t_rep_matrix(raw[[paste0("I", gu)]], n_t, n_rep)
      Iv <- .dt_col_to_t_rep_matrix(raw[[paste0("I", gv)]], n_t, n_rep)
      # Per-capita FOI in site l (matches dust lambda = beta * L * I / sum(N)).
      cum_foi    <- .cum_trapz(beta * (Iu + Iv) / N_site_vec[l], timepoints)
      P_vac_cf   <- 1 - rowMeans(exp(-alpha * cum_foi))
      P_unvac_cf <- 1 - rowMeans(exp(-1     * cum_foi))
      Cu <- rowMeans(.dt_col_to_t_rep_matrix(raw[[paste0("C", gu)]], n_t, n_rep))
      Cv <- rowMeans(.dt_col_to_t_rep_matrix(raw[[paste0("C", gv)]], n_t, n_rep))
      num_t      <- num_t      + Cv + P_vac_cf   * unvac_l[l]
      denom_t    <- denom_t    + Cu + P_unvac_cf * vac_l[l]
      tot_Cvac   <- tot_Cvac   + Cv
      tot_Cunvac <- tot_Cunvac + Cu
    }
    eate_t <- num_t / denom_t
    ave_t  <- (denom_t - num_t) / N_tot
    # Pooled factual CIR across sites.
    P_fac_vac   <- tot_Cvac   / sum(vac_l)
    P_fac_unvac <- tot_Cunvac / sum(unvac_l)
    crr_t     <- P_fac_vac / P_fac_unvac
    crr_ave_t <- P_fac_unvac - P_fac_vac

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = crr_t, ave = crr_ave_t,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  res <- parallel::mclapply(seq_len(n_vac),
                            function(i) run_one_allocation(),
                            mc.cores = mc.cores)
  rbindlist(res, fill = TRUE)
}

# ---------------------------------------------------------------------------
# Two-block "effect modification" model. The population splits into two
# NON-MIXING compartments (block-diagonal contacts): fraction split_frac in A,
# the rest in B. The vaccinated susceptibility differs by compartment --
# alpha in A and split_alpha_prod / alpha in B (so alpha_A * alpha_B =
# split_alpha_prod, fixed) -- i.e. the SAME vaccine has a different effect in
# the two sub-populations. Each realisation is seeded by a single index case
# in a random person, so it lands in A w.p. split_frac (else B) and only that
# compartment ignites. Because A and B have different vaccine ratios, the MODE
# (typical = A outbreak) and the MEAN (A/B blend) give different VE.
#
# Groups: [unvac_A, vac_A, unvac_B, vac_B]; within-compartment mixing weight =
# sum(N)/N_comp so each keeps R0 = beta/gamma. Seeding is stratified (run each
# seed-compartment, weight by split_frac) -- the exact single-seed mixture,
# variance-reduced.
.split_effect_setup <- function(alpha, f, N, split_frac, split_alpha_prod) {
  alpha_B <- split_alpha_prod / alpha
  N_A <- round(split_frac * N); N_B <- N - N_A
  unvacA <- round((1 - f) * N_A); vacA <- N_A - unvacA
  unvacB <- round((1 - f) * N_B); vacB <- N_B - unvacB
  mm <- matrix(0, 4, 4); mm[1:2, 1:2] <- N / N_A; mm[3:4, 3:4] <- N / N_B
  list(Ngrp = c(unvacA, vacA, unvacB, vacB), sus = c(1, alpha, 1, alpha_B),
       mm = mm, N_A = N_A, N_B = N_B,
       unvacA = unvacA, vacA = vacA, unvacB = unvacB, vacB = vacB,
       alpha_A = alpha, alpha_B = alpha_B)
}

get_stoch_eate_sir_split_effect <- function(beta = 1, susceptibility = c(1, 1),
                                            f = 0.5, N = 200, t = 30, gamma = 1,
                                            I_ini_total = 1, split_frac = 0.75,
                                            split_alpha_prod = 0.5,
                                            n_vac = 10, n_rep = 20,
                                            dt = 0.1, timepoints = NULL,
                                            mc.cores = 10, inner_cores = 1,
                                            seed = NULL) {
  alpha <- susceptibility[2]
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t <- length(timepoints)
  s   <- .split_effect_setup(alpha, f, N, split_frac, split_alpha_prod)
  I0  <- as.integer(max(1, round(I_ini_total)))

  # Zero-padded [n_t, n_rep] matrix from a batch that fills columns col_range.
  pad <- function(vec_or_null, n_rep_full, col_range) {
    m <- matrix(0, nrow = n_t, ncol = n_rep_full)
    if (length(col_range)) m[, col_range] <- vec_or_null
    m
  }

  run_one_allocation <- function() {
    sim_id <- runif(1)
    nA <- max(1L, round(split_frac * n_rep)); nB <- n_rep - nA
    run_batch <- function(gseed, nsim) {
      if (nsim <= 0L) return(NULL)
      Ii <- integer(4); Ii[gseed] <- I0
      raw <- run_stoch_cd_dust(s$mm, beta = beta, N = s$Ngrp, t = t, I_ini = Ii,
                               susceptibility = s$sus, gamma = gamma, dt = dt,
                               timepoints = timepoints, n_sim = nsim,
                               cores = inner_cores, seed = seed)
      setDT(raw); raw
    }
    rawA <- run_batch(1L, nA)      # seed in compartment A (group 1)
    rawB <- run_batch(3L, nB)      # seed in compartment B (group 3)

    col <- function(raw, nm, ncol) if (is.null(raw)) NULL else
      .dt_col_to_t_rep_matrix(raw[[nm]], n_t, ncol)
    cA <- if (nB > 0) seq_len(nA) else seq_len(nA)
    cB <- if (nA > 0) (nA + 1L):n_rep else seq_len(nB)

    # Per-compartment infected + case matrices, zero in the non-seeded batch.
    I_A <- pad(if (nA > 0) col(rawA, "I1", nA) + col(rawA, "I2", nA) else NULL, n_rep, cA)
    I_B <- pad(if (nB > 0) col(rawB, "I3", nB) + col(rawB, "I4", nB) else NULL, n_rep, cB)
    Cu_A <- pad(if (nA > 0) col(rawA, "C1", nA) else NULL, n_rep, cA)
    Cv_A <- pad(if (nA > 0) col(rawA, "C2", nA) else NULL, n_rep, cA)
    Cu_B <- pad(if (nB > 0) col(rawB, "C3", nB) else NULL, n_rep, cB)
    Cv_B <- pad(if (nB > 0) col(rawB, "C4", nB) else NULL, n_rep, cB)

    cum_foi_A <- .cum_trapz(beta * I_A / s$N_A, timepoints)
    cum_foi_B <- .cum_trapz(beta * I_B / s$N_B, timepoints)

    # Frozen-field counterfactuals per compartment (each uses its own alpha).
    P_vac_cf_A   <- 1 - rowMeans(exp(-s$alpha_A * cum_foi_A))
    P_unvac_cf_A <- 1 - rowMeans(exp(-1         * cum_foi_A))
    P_vac_cf_B   <- 1 - rowMeans(exp(-s$alpha_B * cum_foi_B))
    P_unvac_cf_B <- 1 - rowMeans(exp(-1         * cum_foi_B))

    P_fac_vac_A   <- rowMeans(Cv_A) / s$vacA
    P_fac_unvac_A <- rowMeans(Cu_A) / s$unvacA
    P_fac_vac_B   <- rowMeans(Cv_B) / s$vacB
    P_fac_unvac_B <- rowMeans(Cu_B) / s$unvacB

    num_t   <- s$vacA   * P_fac_vac_A   + s$unvacA * P_vac_cf_A +
               s$vacB   * P_fac_vac_B   + s$unvacB * P_vac_cf_B
    denom_t <- s$unvacA * P_fac_unvac_A + s$vacA   * P_unvac_cf_A +
               s$unvacB * P_fac_unvac_B + s$vacB   * P_unvac_cf_B
    eate_t  <- num_t / denom_t
    ave_t   <- (denom_t - num_t) / N
    # Pooled factual CIR across compartments.
    P_fac_vac   <- (rowMeans(Cv_A) + rowMeans(Cv_B)) / (s$vacA   + s$vacB)
    P_fac_unvac <- (rowMeans(Cu_A) + rowMeans(Cu_B)) / (s$unvacA + s$unvacB)
    crr_t     <- P_fac_vac / P_fac_unvac
    crr_ave_t <- P_fac_unvac - P_fac_vac

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = crr_t, ave = crr_ave_t,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  res <- parallel::mclapply(seq_len(n_vac),
                            function(i) run_one_allocation(),
                            mc.cores = mc.cores)
  rbindlist(res, fill = TRUE)
}

# Per-node susceptibility for the contact-dependent vaccine effect:
#     alpha_eff[i] = 1 - (f[i] / vac_frac_ref)^vac_frac_power * (1 - alpha)
# with f[i] = (vac[i] + #vaccinated contacts) / (1 + degree[i]), i.e. the
# fraction vaccinated in i's local neighbourhood INCLUDING i itself, and 1 for
# unvaccinated nodes. Counting self is what lets a lone vaccinated person get
# some protection (f = 1/(1+degree) rather than 0). Must match the frac_vac
# definition in stoch_mod_adj_vacfrac.R.
#
# This is a pure function of the parameters -- nothing state-dependent -- so it
# can be precomputed and handed to ANY engine as a plain susceptibility vector.
# That is what lets this model use the event-driven simulator (and makes the
# bespoke odin model optional rather than necessary).
vacfrac_susceptibility <- function(vac, alpha, adj = NULL, c_ij = NULL,
                                   vac_frac_ref = 1, vac_frac_power = 1,
                                   vac_frac_thresh = 0, spill = 0,
                                   trans_tau = 1) {
  if (is.null(adj)) adj <- contact_matrix_to_adj(c_ij)
  n   <- ncol(adj$neighbors)
  vi  <- integer(n); vi[vac] <- 1L
  deg <- colSums(adj$mask)
  nv  <- colSums(adj$mask * matrix(vi[adj$neighbors], nrow(adj$neighbors), n))
  f   <- (vi + nv) / (deg + 1)

  # Protection as a fraction of the full effect (1 - alpha), by local coverage.
  prot <- if (vac_frac_thresh > 0) as.numeric(f >= vac_frac_thresh)
          else (f / vac_frac_ref)^vac_frac_power
  sus <- rep(1, n)
  sus[vac] <- 1 - prot[vac] * (1 - alpha)

  # `spill` (mechanism B): UNVACCINATED people also get protection from living
  # in a well-vaccinated neighbourhood, at fraction `spill` of the vaccinated
  # effect. This is what contaminates a trial's control arm -- a control's
  # outcome now depends on OTHER people's treatment -- so the observed CIR is
  # biased toward 1 while the causal flip-VE is not. spill = 0 recovers the
  # earlier behaviour (unvaccinated susceptibility exactly 1).
  if (spill > 0) {
    unvac <- setdiff(seq_len(n), vac)
    sus[unvac] <- 1 - spill * prot[unvac] * (1 - alpha)
  }

  # `trans_tau` (mechanism A): coverage-dependent TRANSMISSIBILITY. A
  # vaccinated person in a well-vaccinated neighbourhood is also less
  # infectious, trans = 1 - prot * (1 - trans_tau). Flipping i then directly
  # lowers the risk i poses to its contacts -- a spillover that lands inside
  # i's own flip. trans_tau = 1 means no transmissibility effect.
  trans <- rep(1, n)
  if (trans_tau < 1) trans[vac] <- 1 - prot[vac] * (1 - trans_tau)

  attr(sus, "trans") <- trans
  sus
}


# ---------------------------------------------------------------------------
# Network SIR with a CONTACT-DEPENDENT vaccine effect (stoch_mod_adj_vacfrac.R)
# ---------------------------------------------------------------------------
# Same individual-level network model as run_stoch_adj, but a vaccinated node's
# susceptibility depends on the fraction f of ITS OWN contacts that are
# vaccinated:
#     alpha_eff = 1 - f^vac_frac_power * (1 - alpha)
# so f = 0 gives no vaccine effect, f = 1 the full alpha, and vac_frac_power
# bends the interpolation. susceptibility is computed inside the model, so pass
# `vac` (node indices) and the scalar `alpha` rather than a susceptibility
# vector. `adj` may be supplied to skip rebuilding the adjacency list.
run_stoch_adj_vacfrac <- function(contact_matrix, beta, t, I_ini, vac, alpha,
                                  N=NULL, transmissibility=NULL,
                                  vac_frac_power=1, vac_frac_ref=0.5,
                                  gamma=1/3, waning=0, dt=0.1,
                                  timepoints=seq(0, t, 1), n_sim=100, cores=10,
                                  seed=NULL, adj=NULL) {
  adj <- if (is.null(adj)) contact_matrix_to_adj(contact_matrix) else adj
  n   <- if (is.null(contact_matrix)) ncol(adj$neighbors) else nrow(contact_matrix)

  if (is.null(N))                N                <- rep(1L, n)
  if (is.null(transmissibility)) transmissibility <- rep(1, n)
  vac_ind      <- integer(n); vac_ind[vac] <- 1L

  beta_scalar <- if (length(beta) == 1L) beta else mean(beta)

  params <- list(
    n               = as.integer(n),
    max_degree      = as.integer(adj$max_degree),
    beta            = beta_scalar,
    gamma           = gamma,
    waning          = waning,
    alpha           = alpha,
    vac_frac_power  = vac_frac_power,
    vac_frac_ref    = vac_frac_ref,
    neighbors       = adj$neighbors,
    mask            = adj$mask,
    vac             = vac_ind,
    transmisibility = transmissibility,
    S_ini           = N - I_ini,
    I_ini           = I_ini
  )

  sys <- dust2::dust_system_create(stoch_model_adj_vacfrac, params,
                                   n_particles = n_sim,
                                   n_threads   = cores,
                                   dt          = dt,
                                   time        = 0,
                                   seed        = seed)
  dust2::dust_system_set_state_initial(sys)
  raw <- dust2::dust_system_simulate(sys, timepoints)
  if (length(dim(raw)) == 2L) {
    raw <- array(raw, dim = c(dim(raw)[1L], 1L, dim(raw)[2L]))
  }
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


# run_stoch_network counterpart: returns (time, sim, C1 unvaccinated cases,
# C2 vaccinated cases). beta is on the user scale (multiplied by N/k_mean
# internally, as in run_stoch_network).
run_stoch_network_vacfrac <- function(beta=1, N=100, pl_alpha=3, alpha=0.5,
                                      t=100, vac_frac=0.5, vac=NULL,
                                      vac_frac_power=1, vac_frac_ref=0.5,
                                      gamma=1/3,
                                      c_ij=NULL, k_mean=6,
                                      dt=0.1, timepoints=seq(0, t, 1), I_ini=2,
                                      n_sim=100, cores=10, seed=NULL, adj=NULL,
                                      engine=c("dust", "events"), csr=NULL) {
  engine <- match.arg(engine)
  if (is.null(c_ij) && is.null(adj) && is.null(csr))
    c_ij <- get_conact_matrix_pl(N, pl_alpha, mean_k=k_mean)
  if (is.null(vac)) vac <- sample(seq_len(N), vac_frac * N)
  non_vac <- setdiff(seq_len(N), vac)

  # alpha_eff is a precomputable per-node vector, so the exact event-driven
  # engine can run this model unchanged -- ~137x faster than tau-leaping, and
  # without dust's dt bias.
  if (engine == "events") {
    if (is.null(csr)) csr <- adj_to_csr(contact_matrix = c_ij, adj = adj)
    sus <- vacfrac_susceptibility(vac, alpha, adj = adj, c_ij = c_ij,
                                  vac_frac_ref = vac_frac_ref,
                                  vac_frac_power = vac_frac_power)
    return(run_stoch_network_events(
      beta = beta, N = N, susceptibility = sus, t = t, vac = vac, csr = csr,
      gamma = gamma, timepoints = timepoints, I_ini = I_ini, n_sim = n_sim,
      seed = seed, k_mean = k_mean, cores = cores))
  }

  I_ini <- c(rep(1, I_ini), rep(0L, N - I_ini))

  full <- run_stoch_adj_vacfrac(c_ij, beta = N * beta / k_mean, t = t,
                                I_ini = I_ini, vac = vac, alpha = alpha,
                                vac_frac_power = vac_frac_power,
                                vac_frac_ref = vac_frac_ref,
                                gamma = gamma, dt = dt, timepoints = timepoints,
                                n_sim = n_sim, cores = cores, seed = seed,
                                adj = adj)

  vac_cols   <- paste0("C", vac)
  unvac_cols <- paste0("C", non_vac)
  full[, .(time, sim,
           C2 = rowSums(.SD[, vac_cols,   with=FALSE]),
           C1 = rowSums(.SD[, unvac_cols, with=FALSE]))]
}

# ---------------------------------------------------------------------------
# Interference-aware EATE for the contact-dependent-vaccine-effect network
# ---------------------------------------------------------------------------
# get_stoch_eate_network uses a FROZEN FIELD: the counterfactual is obtained by
# re-weighting the realised force of infection, which assumes flipping one
# person does not change the epidemic. That fails here -- flipping i changes
# f[j] for every neighbour j (by 1/degree(j), ~0.17 at mean degree 6), hence
# their protection, hence the FOI. So the counterfactual arm is obtained by
# RE-SIMULATING the model with i's status flipped, which captures that
# feedback (including back onto i).
#
# Only the flipped arm needs simulating: the factual run already supplies each
# individual's observed outcome. Cost is therefore 1 + n_flip runs per
# allocation, not 2 * n_flip.
#
# All runs share a seed (common random numbers), so the factual and flipped
# runs are coupled and the difference is not swamped by simulation noise.
#
# n_flip individuals are sampled per allocation (stratified across the
# vaccinated and unvaccinated arms) and their group means scaled to the group
# sizes, so n_flip trades accuracy against cost; n_flip = N is exact.
get_stoch_eate_network_vacfrac <- function(beta = 1, alpha = 0.5, f = 0.5,
                                           N = 200, t = 15, pl_alpha = 3,
                                           c_ij = NULL, vac_frac_power = 1,
                                           vac_frac_ref = 1, vac_frac_thresh = 0,
                                           n_vac = 10, n_rep = 20, n_flip = 20,
                                           k_mean = 6, gamma = 1 / 3, dt = 0.1,
                                           timepoints = NULL, init_I = 2,
                                           mc.cores = 10, inner_cores = 1,
                                           vac_list = NULL, adj = NULL,
                                           crn_seed = 1L,
                                           engine = c("dust", "events"),
                                           csr = NULL) {
  engine <- match.arg(engine)
  if (is.null(c_ij) && is.null(adj))
    c_ij <- get_conact_matrix_pl(N, pl_alpha, mean_k = k_mean)
  if (is.null(adj)) adj <- contact_matrix_to_adj(c_ij)
  if (engine == "events" && is.null(csr))
    csr <- adj_to_csr(contact_matrix = c_ij, adj = adj)
  if (is.null(timepoints)) timepoints <- seq(1, t, 1)
  n_t <- length(timepoints)
  if (!is.null(vac_list)) n_vac <- length(vac_list)
  I_ini_vec <- c(rep(1L, init_I), rep(0L, N - init_I))

  run_one_allocation <- function(a) {
    vac     <- if (!is.null(vac_list)) vac_list[[a]]
               else sample(seq_len(N), round(f * N))
    non_vac <- setdiff(seq_len(N), vac)
    sim_id  <- runif(1)
    seed_a  <- as.integer(crn_seed) + a          # CRN: shared by all runs here

    # Both engines return a [n_t, N] matrix of per-individual P(infected by t),
    # which is what the counterfactual contrast needs.
    sim_v <- if (engine == "events") function(v) {
      sus <- vacfrac_susceptibility(v, alpha, adj = adj,
                                    vac_frac_ref = vac_frac_ref,
                                    vac_frac_power = vac_frac_power,
                                    vac_frac_thresh = vac_frac_thresh)
      inf <- run_stoch_network_events(
        beta = beta, N = N, susceptibility = sus, t = t, vac = v, csr = csr,
        gamma = gamma, timepoints = timepoints, I_ini = init_I, n_sim = n_rep,
        seed = seed_a, k_mean = k_mean, cores = inner_cores,
        return_times = TRUE)                     # [n_rep, N] infection times
      vapply(timepoints, function(tp) colMeans(inf <= tp), numeric(N)) |> t()
    } else function(v) {
      r <- run_stoch_adj_vacfrac(c_ij, beta = N * beta / k_mean, t = t,
                                 I_ini = I_ini_vec, vac = v, alpha = alpha,
                                 vac_frac_power = vac_frac_power,
                                 vac_frac_ref = vac_frac_ref,
                                 gamma = gamma, dt = dt, timepoints = timepoints,
                                 n_sim = n_rep, cores = inner_cores,
                                 seed = seed_a, adj = adj)
      setDT(r)
      vapply(seq_len(N), function(k)
        rowMeans(.dt_col_to_t_rep_matrix(r[[paste0("C", k)]], n_t, n_rep)),
        numeric(n_t))
    }

    P_fac <- sim_v(vac)                          # [n_t, N]

    # Sample individuals to flip, stratified across the two arms.
    n_fv   <- max(1L, round(n_flip * length(vac) / N))
    n_fu   <- max(1L, n_flip - n_fv)
    flip_v <- if (length(vac))     sample(vac,     min(n_fv, length(vac)))     else integer(0)
    flip_u <- if (length(non_vac)) sample(non_vac, min(n_fu, length(non_vac))) else integer(0)

    # Counterfactual by RE-SIMULATION: flip i, take that run's mean for i.
    cf <- function(i, vaccinate) {
      v2 <- if (vaccinate) sort(c(vac, i)) else setdiff(vac, i)
      sim_v(v2)[, i]
    }
    P_cf_u <- if (length(flip_u))                 # unvaccinated, if vaccinated
      vapply(flip_u, function(i) cf(i, TRUE),  numeric(n_t)) else NULL
    P_cf_v <- if (length(flip_v))                 # vaccinated, if unvaccinated
      vapply(flip_v, function(i) cf(i, FALSE), numeric(n_t)) else NULL

    mrow <- function(m) if (is.null(m)) rep(NA_real_, n_t) else
      if (is.matrix(m)) rowMeans(m) else m

    # Same hybrid contrast as get_stoch_eate_network: factual for the matching
    # arm, counterfactual for the flipped arm -- but re-simulated, and the
    # sampled group means scaled up to the group sizes.
    num_t   <- rowSums(P_fac[, vac,     drop = FALSE]) + length(non_vac) * mrow(P_cf_u)
    denom_t <- rowSums(P_fac[, non_vac, drop = FALSE]) + length(vac)     * mrow(P_cf_v)
    eate_t  <- num_t / denom_t
    ave_t   <- (denom_t - num_t) / N

    ar_fac_vac   <- rowSums(P_fac[, vac,     drop = FALSE]) / max(length(vac), 1)
    ar_fac_unvac <- rowSums(P_fac[, non_vac, drop = FALSE]) / max(length(non_vac), 1)

    rbindlist(list(
      data.frame(t = timepoints, eate = eate_t, ave = ave_t,
                 num = num_t, denom = denom_t,
                 method = "full_stoch", sim = sim_id),
      data.frame(t = timepoints, eate = ar_fac_vac / ar_fac_unvac,
                 ave = ar_fac_unvac - ar_fac_vac,
                 num = NA_real_, denom = NA_real_,
                 method = "CRR", sim = sim_id)
    ), fill = TRUE)
  }

  rbindlist(parallel::mclapply(seq_len(n_vac), run_one_allocation,
                               mc.cores = mc.cores), fill = TRUE)
}
