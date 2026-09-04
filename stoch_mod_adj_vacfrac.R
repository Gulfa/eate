# Stochastic adjacency-list SIR with CONTACT-DEPENDENT vaccine effect.
#
# Identical to stoch_mod_adj.R except that a vaccinated node's susceptibility
# is not a fixed alpha: it depends on how many of that node's CONTACTS are
# vaccinated. Writing f[i] for the fraction of i's contacts that are
# vaccinated,
#
#     alpha_eff[i] = alpha ^ ((f[i] / vac_frac_ref) ^ vac_frac_power)
#
# normalised so the model AGREES WITH THE PLAIN NETWORK MODEL at the reference
# coverage (default 0.5, the coverage these runs use):
#     f = 0            -> alpha_eff = 1       (no vaccine effect at all)
#     f = vac_frac_ref -> alpha_eff = alpha   (as in stoch_mod_adj.R)
#     f = 1 (ref 0.5)  -> alpha_eff = alpha^2
#
# Unvaccinated nodes always have susceptibility 1.
#
# Unlike stoch_mod_adj.R, susceptibility is COMPUTED here rather than passed
# in, so the caller supplies vac[] (0/1 per node) and the scalar alpha.

n          <- parameter()
max_degree <- parameter()
beta       <- parameter()
gamma      <- parameter(0.1)
waning     <- parameter(0)

alpha           <- parameter()
vac_frac_power  <- parameter(1)
vac_frac_ref    <- parameter(0.5)   # coverage at which alpha_eff == alpha

dim(neighbors) <- c(max_degree, n)
neighbors <- parameter(type = "integer")

dim(mask)            <- c(max_degree, n)
dim(vac)             <- n
dim(transmisibility) <- n
dim(S_ini)           <- n
dim(I_ini)           <- n
mask            <- parameter()
vac             <- parameter(type = "integer")
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

# --- contact-dependent vaccine effect (static: derived from parameters) -----
degree[] <- sum(mask[, i])
dim(degree) <- n

vac_contrib[, ] <- mask[i, j] * vac[neighbors[i, j]]
dim(vac_contrib) <- c(max_degree, n)

# Fraction of i's contacts that are vaccinated (0 for isolated nodes).
frac_vac[] <- if (degree[i] > 0) sum(vac_contrib[, i]) / degree[i] else 0
dim(frac_vac) <- n

# Vaccinated nodes get the coverage-scaled effect; unvaccinated get none.
#
#     alpha_eff[i] = alpha ^ ((f[i] / vac_frac_ref) ^ vac_frac_power)
#
# Normalised so that AT THE REFERENCE COVERAGE the model reproduces the plain
# network model: a node whose contacts are vac_frac_ref vaccinated gets exactly
# alpha, and with vac_frac_power = 1 the exponent is linear in f, so the
# geometric mean of alpha_eff over vaccinated nodes is exactly alpha whenever
# mean(f) = vac_frac_ref (which is the case at that coverage). Set
# vac_frac_ref = 0.5 to match the usual 50% coverage runs.
#
#   f = 0            -> alpha^0 = 1        (no vaccine effect)
#   f = vac_frac_ref -> alpha              (the standard model's effect)
#   f = 1 (ref 0.5)  -> alpha^2            (stronger where contacts are all vaccinated)
#
# Multiplicative rather than the linear 1 - (f/ref)(1-alpha): the linear form
# goes NEGATIVE once f > ref/(1-alpha) (at ref=0.5, alpha=0.3 that is f > 0.71,
# which happens for a few percent of nodes at mean degree 6), and clamping
# there would distort the very average we are trying to fix.
susceptibility[] <- if (vac[i] == 1)
                      alpha^((frac_vac[i] / vac_frac_ref)^vac_frac_power) else 1
dim(susceptibility) <- n

# --- transmission (unchanged from stoch_mod_adj.R) --------------------------
contrib[, ] <- mask[i, j] * I[neighbors[i, j]] * transmisibility[neighbors[i, j]]
dim(contrib) <- c(max_degree, n)

foi[] <- beta / N_total * sum(contrib[, i])
dim(foi) <- n

N_total <- sum(S) + sum(I) + sum(R)

p_SI[] <- 1 - exp(-susceptibility[i] * foi[i] * dt)
p_IR   <- 1 - exp(-gamma  * dt)
p_RS   <- 1 - exp(-waning * dt)
dim(p_SI) <- n

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
