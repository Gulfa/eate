# Stochastic adjacency-list SIR with CONTACT-DEPENDENT vaccine effect.
#
# Identical to stoch_mod_adj.R except that a vaccinated node's susceptibility
# is not a fixed alpha: it depends on how many of that node's CONTACTS are
# vaccinated. Writing f[i] for the fraction vaccinated in i's local
# neighbourhood -- i itself plus its contacts --
#
#     alpha_eff[i] = 1 - (f[i] / vac_frac_ref) ^ vac_frac_power * (1 - alpha)
#
# i.e. the vaccine EFFECT (1 - alpha) is scaled by the local coverage, reaching
# its full value when the neighbourhood is fully vaccinated:
#     f -> 0            -> alpha_eff -> 1        (no vaccine effect at all)
#     f = 1/(1+degree)  -> 1/(1+degree) of the effect (a lone vaccinated person)
#     f = 1             -> alpha_eff = alpha     (the full effect)
# With vac_frac_ref = 1 (the default) and f in [0, 1], alpha_eff stays in
# [alpha, 1]: it can never overshoot the full effect nor go negative. Setting
# vac_frac_ref < 1 reaches the full effect earlier, but then f > vac_frac_ref
# would overshoot, so keep it at 1 unless you mean that.
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
vac_frac_ref    <- parameter(1)     # local coverage at which the FULL effect is reached
vac_frac_thresh <- parameter(0)     # >0: hard threshold at this local coverage

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

# Fraction vaccinated in i's LOCAL NEIGHBOURHOOD, counting i itself:
#     f[i] = (vac[i] + #vaccinated contacts) / (1 + degree[i])
# Including self is what gives a lone vaccinated person some protection: with
# no vaccinated contacts f = 1/(1+degree) rather than 0, so alpha_eff < 1. It
# also means f[i] depends on i's OWN status, so flipping i changes i's own
# effect -- the local effect a sparse trial should still see. The denominator
# is >= 1, so isolated nodes need no special case.
frac_vac[] <- (vac[i] + sum(vac_contrib[, i])) / (degree[i] + 1)
dim(frac_vac) <- n

# Vaccinated nodes get the coverage-scaled effect; unvaccinated get none.
#
#     alpha_eff[i] = 1 - (f[i] / vac_frac_ref)^vac_frac_power * (1 - alpha)
#
# The vaccine EFFECT (1 - alpha) is scaled by local coverage: a lone vaccinated
# person (f = 1/(1+degree)) gets 1/(1+degree) of the effect, a fully vaccinated
# neighbourhood (f = 1) gets all of it. Linear in f with vac_frac_ref = 1 and
# f in [0, 1] keeps alpha_eff in [alpha, 1] -- it cannot overshoot the full
# effect, and cannot go negative (which is what ruled the linear form out when
# the reference was 0.5 rather than 1).
# vac_frac_thresh > 0 switches from the smooth ramp to a HARD THRESHOLD:
# full effect (alpha) once local coverage reaches vac_frac_thresh, none below.
# A step makes a single flip able to push a neighbour ACROSS the threshold, so
# one vaccination can discontinuously change someone else's protection.
susceptibility[] <- if (vac[i] == 0) 1 else if (vac_frac_thresh > 0)
                     (if (frac_vac[i] >= vac_frac_thresh) alpha else 1) else
                     1 - (frac_vac[i] / vac_frac_ref)^vac_frac_power * (1 - alpha)
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
