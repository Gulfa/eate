# Stochastic SIR model with n compartments — discrete-time tau-leaping.
# odin2 / dust2 model. Parameter interface mirrors det_mod_cd.R (same names,
# same shapes), plus a step size `dt`. Replication is handled by dust2 via
# n_particles, so this model file does not loop over simulations itself.

n       <- parameter()
N_steps <- parameter()
gamma   <- parameter(0.1)
waning  <- parameter(0)

dim(beta_day) <- N_steps
beta_day      <- parameter()
# One beta value per dt step. At simulated time `time` (start of step) we read
# beta_day[floor(time/dt) + 1]; first step (time=0) reads beta_day[1].
beta <- beta_day[as.integer(time / dt) + 1]

dim(mixing_matrix)   <- c(n, n)
dim(S_ini)           <- n
dim(I_ini)           <- n
dim(susceptibility)  <- n
dim(transmisibility) <- n
mixing_matrix   <- parameter()
S_ini           <- parameter()
I_ini           <- parameter()
susceptibility  <- parameter()
transmisibility <- parameter()

# State
initial(S[])        <- S_ini[i]
initial(I[])        <- I_ini[i]
initial(R[])        <- 0
initial(C[])        <- 0
initial(exposure[]) <- 0

dim(S)        <- n
dim(I)        <- n
dim(R)        <- n
dim(C)        <- n
dim(N)        <- n
dim(exposure) <- n

# Force of infection (group i acquires from contacts with group j)
lambda_ij[, ] <- beta * mixing_matrix[i, j] * I[j] / sum(N) * transmisibility[j]
dim(lambda_ij) <- c(n, n)

# Per-step transition probabilities
p_SI[] <- 1 - exp(-susceptibility[i] * sum(lambda_ij[i, ]) * dt)
p_IR   <- 1 - exp(-gamma  * dt)
p_RS   <- 1 - exp(-waning * dt)
dim(p_SI) <- n

# Binomial draws
n_SI[]     <- Binomial(S[i], p_SI[i])
n_IR[]     <- Binomial(I[i], p_IR)
n_RS[]     <- Binomial(R[i], p_RS)
n_RS_in[]  <- 0
n_RS_in[1] <- sum(n_RS[])
dim(n_SI)    <- n
dim(n_IR)    <- n
dim(n_RS)    <- n
dim(n_RS_in) <- n

# Compartment updates
update(S[])        <- S[i] - n_SI[i] + n_RS_in[i]
update(I[])        <- I[i] + n_SI[i] - n_IR[i]
update(R[])        <- R[i] + n_IR[i] - n_RS[i]
update(C[])        <- C[i] + n_SI[i]
update(exposure[]) <- exposure[i] + sum(lambda_ij[i, ]) * dt

# Per-group total population
N[] <- S[i] + I[i] + R[i]
