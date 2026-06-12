# Stochastic linear (exposure-only) model with n independent groups.
# Each S[i] becomes C with rate susceptibility[i] * beta per individual per
# unit time — no interaction between groups, no recovery. Discrete-time
# tau-leaping with binomial draws over step `dt` (provided by dust2).

n    <- parameter()
beta <- parameter()

dim(S_ini)          <- n
dim(susceptibility) <- n
S_ini          <- parameter()
susceptibility <- parameter()

initial(S[]) <- S_ini[i]
initial(C[]) <- 0
dim(S) <- n
dim(C) <- n

# Per-step transition probability
p_SC[] <- 1 - exp(-susceptibility[i] * beta * dt)
dim(p_SC) <- n

# Binomial draw for new exposures
n_SC[] <- Binomial(S[i], p_SC[i])
dim(n_SC) <- n

update(S[]) <- S[i] - n_SC[i]
update(C[]) <- C[i] + n_SC[i]
