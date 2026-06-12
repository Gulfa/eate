# Individual-level stochastic SIR — discrete-time tau-leaping (odin2/dust2).
#
# Each person i (1..n) carries integer S[i], I[i], R[i] in {0, 1} (with
# S[i] + I[i] + R[i] = 1). Contacts are stored as a weighted sparse
# adjacency:
#   neighbors[k, i] = identity of person i's k-th contact (padded with 1)
#   weights[k, i]   = strength of that contact (0 for padding slots)
# The FOI on person i is
#   beta * sum_k weights[k, i] * I[neighbors[k, i]] * transmisibility[neighbors[k, i]] / N_total
# `foi[i]` is stored as a state so it is recorded at every requested time;
# foi[i] at time t holds the FOI that drove person i's transition during
# the step ending at t (i.e. computed from I at time t - dt).

n          <- parameter()
max_degree <- parameter()
beta       <- parameter()
gamma      <- parameter(0.1)

dim(neighbors) <- c(max_degree, n)
neighbors      <- parameter(type = "integer")
dim(weights)   <- c(max_degree, n)
weights        <- parameter()

dim(susceptibility)  <- n
dim(transmisibility) <- n
dim(S_ini)           <- n
dim(I_ini)           <- n
susceptibility  <- parameter()
transmisibility <- parameter()
S_ini           <- parameter()
I_ini           <- parameter()

initial(S[])   <- S_ini[i]
initial(I[])   <- I_ini[i]
initial(R[])   <- 0
initial(foi[]) <- 0
dim(S)   <- n
dim(I)   <- n
dim(R)   <- n
dim(foi) <- n

# Per-slot contribution to person i's force of infection
contrib[, ] <- weights[i, j] * I[neighbors[i, j]] * transmisibility[neighbors[i, j]]
dim(contrib) <- c(max_degree, n)

# Current-step FOI per person (column sum — requires odin2)
foi_now[] <- beta * sum(contrib[, i]) / N_total
dim(foi_now) <- n

N_total <- sum(S) + sum(I) + sum(R)

# Transition probabilities over a step of length dt
p_SI[] <- 1 - exp(-susceptibility[i] * foi_now[i] * dt)
p_IR   <- 1 - exp(-gamma * dt)
dim(p_SI) <- n

# Binomial draws — S[i] / I[i] are 0 or 1, so these are Bernoulli per person.
n_SI[] <- Binomial(S[i], p_SI[i])
n_IR[] <- Binomial(I[i], p_IR)
dim(n_SI) <- n
dim(n_IR) <- n

update(S[])   <- S[i] - n_SI[i]
update(I[])   <- I[i] + n_SI[i] - n_IR[i]
update(R[])   <- R[i] + n_IR[i]
update(foi[]) <- foi[i] + foi_now[i]
