# Tests for the event-driven network SIR (net_sir_events.cpp) and for the
# engine switch that selects it inside run_stoch_network().

skip_if_not_installed("data.table")
skip_if_not_installed("cpp11")

library(data.table)

if (!exists("run_stoch_network")) {
  source("../../utils.R", chdir=TRUE)
  source("../../stoch_model.R", chdir=TRUE)
}

# ---------------------------------------------------------------------------
# Engine selection
# ---------------------------------------------------------------------------

test_that("network_engine() honours the option and the env var", {
  old_opt <- getOption("eate.network_engine")
  old_env <- Sys.getenv("EATE_NETWORK_ENGINE", unset = NA)
  on.exit({
    options(eate.network_engine = old_opt)
    if (is.na(old_env)) Sys.unsetenv("EATE_NETWORK_ENGINE")
    else Sys.setenv(EATE_NETWORK_ENGINE = old_env)
  }, add = TRUE)

  Sys.unsetenv("EATE_NETWORK_ENGINE")
  options(eate.network_engine = NULL)
  expect_equal(network_engine(), "events")          # the new default

  options(eate.network_engine = "dust")
  expect_equal(network_engine(), "dust")

  Sys.setenv(EATE_NETWORK_ENGINE = "events")        # env var wins over option
  expect_equal(network_engine(), "events")

  Sys.setenv(EATE_NETWORK_ENGINE = "nonsense")
  expect_error(network_engine())
})

test_that("both network engines return the same columns and shape", {
  set.seed(11)
  N <- 60
  cm  <- get_conact_matrix_pl(N, alpha=3, mean_k=4)
  adj <- contact_matrix_to_adj(cm)
  vac <- sample(seq_len(N), N/2)
  tp  <- seq(1, 5, 1)
  args <- list(beta=2, N=N, susceptibility=c(1, 0.5), t=5, c_ij=cm, adj=adj,
               vac=vac, k_mean=4, gamma=1, dt=0.05, timepoints=tp, I_ini=2,
               n_sim=8, cores=1, seed=3)
  d <- do.call(run_stoch_network, c(args, list(engine="dust")))
  e <- do.call(run_stoch_network, c(args, list(engine="events")))
  expect_identical(names(d), c("time", "sim", "C1", "C2"))
  expect_identical(names(e), names(d))
  expect_identical(dim(d), dim(e))
  expect_setequal(e$time, tp)
  expect_setequal(e$sim, seq_len(8))
})

# ---------------------------------------------------------------------------
# The two engines must count the same thing
# ---------------------------------------------------------------------------

test_that("neither network engine counts the index cases (beta = 0)", {
  # stoch_mod_adj.R has initial(C[]) <- 0, so the I_ini seeds are never in C.
  # The event simulator records an infection time of 0 for each seed, so its
  # wrapper has to drop them; otherwise every count carries a constant +I_ini.
  set.seed(12)
  N <- 60
  cm  <- get_conact_matrix_pl(N, alpha=3, mean_k=4)
  adj <- contact_matrix_to_adj(cm)
  vac <- sample(seq_len(N), N/2)
  args <- list(beta=0, N=N, susceptibility=c(1, 1), t=5, c_ij=cm, adj=adj,
               vac=vac, k_mean=4, gamma=1, dt=0.05, timepoints=seq(1, 5, 1),
               I_ini=3, n_sim=8, cores=1, seed=3)
  d <- do.call(run_stoch_network, c(args, list(engine="dust")))
  e <- do.call(run_stoch_network, c(args, list(engine="events")))
  expect_true(all(d$C1 + d$C2 == 0))
  expect_true(all(e$C1 + e$C2 == 0))

  # ... and counting them is still available explicitly, for anyone who wants
  # incidence including the seeds.
  csr <- adj_to_csr(adj = adj)
  s <- run_stoch_network_events(beta=0, N=N, susceptibility=c(1, 1), t=5,
                                vac=vac, csr=csr, gamma=1,
                                timepoints=seq(1, 5, 1), I_ini=3, n_sim=8,
                                seed=3, k_mean=4, cores=1, count_seeds=TRUE)
  expect_true(all(s$C1 + s$C2 == 3))
})

test_that("event engine agrees with dust on the mean at small dt", {
  # Same continuous-time process: the tau-leaping bias vanishes as dt -> 0, so
  # the two must agree to within Monte-Carlo error on a small graph.
  set.seed(15)
  N <- 120
  cm  <- get_conact_matrix_pl(N, alpha=3, mean_k=5)
  adj <- contact_matrix_to_adj(cm)
  vac <- sample(seq_len(N), N/2)
  tp  <- seq(1, 6, 1)
  args <- list(beta=2.0, N=N, susceptibility=c(1, 0.4), t=6, c_ij=cm, adj=adj,
               vac=vac, k_mean=5, gamma=1, timepoints=tp, I_ini=2,
               n_sim=600, cores=1, seed=21)
  d <- do.call(run_stoch_network, c(args, list(engine="dust", dt=0.005)))
  e <- do.call(run_stoch_network, c(args, list(engine="events")))
  fd <- d[time == 6]; fe <- e[time == 6]
  z <- function(x, y) abs(mean(x) - mean(y)) /
    sqrt(var(x)/length(x) + var(y)/length(y))
  expect_lt(z(fd$C1, fe$C1), 4)
  expect_lt(z(fd$C2, fe$C2), 4)
})

# ---------------------------------------------------------------------------
# Determinism and graph construction
# ---------------------------------------------------------------------------

test_that("event engine is reproducible and thread-count independent", {
  set.seed(13)
  N <- 80
  cm  <- get_conact_matrix_pl(N, alpha=3, mean_k=5)
  csr <- adj_to_csr(contact_matrix = cm)
  vac <- sample(seq_len(N), N/2)
  args <- list(beta=2.5, N=N, susceptibility=c(1, 0.4), t=6, vac=vac, csr=csr,
               gamma=1, timepoints=seq(1, 6, 1), I_ini=2, n_sim=16, seed=99,
               k_mean=5)
  a  <- do.call(run_stoch_network_events, c(args, list(cores=1)))
  b  <- do.call(run_stoch_network_events, c(args, list(cores=1)))
  c4 <- do.call(run_stoch_network_events, c(args, list(cores=4)))
  expect_equal(a, b)
  expect_equal(a, c4)          # per-realisation RNG streams, not per-thread
})

test_that("an unset seed still moves with the ambient RNG", {
  set.seed(16)
  N <- 80
  cm  <- get_conact_matrix_pl(N, alpha=3, mean_k=5)
  csr <- adj_to_csr(contact_matrix = cm)
  vac <- sample(seq_len(N), N/2)
  run <- function() run_stoch_network_events(
    beta=2.5, N=N, susceptibility=c(1, 0.4), t=6, vac=vac, csr=csr, gamma=1,
    timepoints=seq(1, 6, 1), I_ini=2, n_sim=16, seed=NULL, k_mean=5, cores=1)
  set.seed(77); a <- run()
  set.seed(77); b <- run()
  d <- run()
  expect_equal(a, b)                       # set.seed() controls it
  expect_false(isTRUE(all.equal(a, d)))     # consecutive calls differ
})

test_that("event engine can build its own graph from adj or c_ij", {
  set.seed(14)
  N <- 60
  cm  <- get_conact_matrix_pl(N, alpha=3, mean_k=4)
  adj <- contact_matrix_to_adj(cm)
  vac <- sample(seq_len(N), N/2)
  args <- list(beta=2, N=N, susceptibility=c(1, 0.5), t=5, vac=vac, gamma=1,
               timepoints=seq(1, 5, 1), I_ini=2, n_sim=8, seed=5, k_mean=4,
               cores=1)
  from_csr <- do.call(run_stoch_network_events,
                      c(args, list(csr=adj_to_csr(adj = adj))))
  from_adj <- do.call(run_stoch_network_events, c(args, list(adj=adj)))
  from_cij <- do.call(run_stoch_network_events, c(args, list(c_ij=cm)))
  expect_equal(from_csr, from_adj)
  expect_equal(from_csr, from_cij)
  expect_error(do.call(run_stoch_network_events, args))   # no graph at all
})

test_that("adj_to_csr transposes into the transmit direction", {
  # get_conact_matrix_pl is directed: row i lists who infects i. The event
  # simulator walks outward from an infectious node, so it needs column k.
  cm <- matrix(0L, 3, 3)
  cm[2, 1] <- 1L        # node 1 infects node 2
  cm[3, 1] <- 1L        # node 1 infects node 3
  csr <- adj_to_csr(contact_matrix = cm)
  out_of <- function(k) sort(csr$nbr[(csr$ptr[k] + 1L):csr$ptr[k + 1L]]) + 1L
  expect_equal(csr$ptr[2] - csr$ptr[1], 2L)   # node 1 has 2 outgoing edges
  expect_equal(out_of(1), c(2L, 3L))
  expect_equal(csr$ptr[3] - csr$ptr[2], 0L)   # node 2 transmits to nobody
  expect_equal(csr$ptr[4] - csr$ptr[3], 0L)
})

test_that("an isolated node is never infected by the event engine", {
  N <- 40
  cm <- matrix(0L, N, N)
  for (i in 2:(N - 1)) { cm[i, i - 1] <- 1L; cm[i - 1, i] <- 1L }  # node N isolated
  csr <- adj_to_csr(contact_matrix = cm)
  ev <- run_stoch_network_events(beta=20, N=N, susceptibility=c(1, 1), t=20,
                                 vac=integer(0), csr=csr, gamma=0.2,
                                 timepoints=20, I_ini=1, n_sim=20, seed=4,
                                 k_mean=2, cores=1)
  expect_true(all(ev$C1 <= N - 2))   # never the seed, never the isolated node
})
