# Stochastic adjacency-list SIR model — discrete-time tau-leaping.
# odin2 / dust2 counterpart of det_mod_adj.R. Same sparse FOI computation
# (O(max_degree) per node) via precomputed neighbors[max_degree, n] and
# mask[max_degree, n]; same parameter interface; deriv() replaced by
# update() with binomial draws over a step of size dt.

n          <- parameter()
max_degree <- parameter()
beta       <- parameter()
gamma      <- parameter(0.1)
waning     <- parameter(0)

dim(neighbors) <- c(max_degree, n)
neighbors <- parameter(type = "integer")

dim(mask)            <- c(max_degree, n)
dim(susceptibility)  <- n
dim(transmisibility) <- n
dim(S_ini)           <- n
dim(I_ini)           <- n
mask            <- parameter()
susceptibility  <- parameter()
transmisibility <- parameter()
S_ini           <- parameter()
I_ini           <- parameter()

initial(S[]) <- S_ini[i]
initial(I[]) <- I_ini[i]
initial(R[]) <- 0
initial(C[]) <- 0
dim(S) <- n
dim(I) <- n
dim(R) <- n
dim(C) <- n

# Per-slot infectivity contribution (max_degree x n)
contrib[, ] <- mask[i, j] * I[neighbors[i, j]] * transmisibility[neighbors[i, j]]
dim(contrib) <- c(max_degree, n)

# FOI per node — column sum (sum over neighbour slots for fixed node i).
# This is the odin2-specific construct that motivated det_mod_adj.R.
foi[] <- beta / N_total * sum(contrib[, i])
dim(foi) <- n

N_total <- sum(S) + sum(I) + sum(R)

# Per-step transition probabilities
p_SI[] <- 1 - exp(-susceptibility[i] * foi[i] * dt)
p_IR   <- 1 - exp(-gamma  * dt)
p_RS   <- 1 - exp(-waning * dt)
dim(p_SI) <- n

# Binomial draws — works for both individual-level (S[i] in {0,1}) and
# grouped (S[i] >= 1) parameterisations.
n_SI[] <- Binomial(S[i], p_SI[i])
n_IR[] <- Binomial(I[i], p_IR)
n_RS[] <- Binomial(R[i], p_RS)
dim(n_SI) <- n
dim(n_IR) <- n
dim(n_RS) <- n

update(S[]) <- S[i] - n_SI[i] + n_RS[i]
update(I[]) <- I[i] + n_SI[i] - n_IR[i]
update(R[]) <- R[i] + n_IR[i] - n_RS[i]
update(C[]) <- C[i] + n_SI[i]
