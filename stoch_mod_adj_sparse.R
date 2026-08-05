# Stochastic adjacency-list SIR model with GROUP-AGGREGATE output. Dynamics
# are identical to stoch_mod_adj.R (per-node discrete-time tau-leaping over
# the adjacency structure). The difference is that the model additionally
# maintains n_groups scalar aggregates (S_g, I_g, R_g, C_g) and dust2 is
# asked to return only those, so the output buffer scales with n_groups
# instead of N. This is the Covasim-style pattern: full per-agent state
# lives inside dust, but only the summaries requested at construction time
# are exposed.
#
# Groups are user-defined via a group indicator matrix G[g, i] = 1 if node
# i belongs to group g, else 0. Groups can overlap (e.g. vac ⊂ pop, ctrl ⊂
# unvac) or partition — the model just sums.

n          <- parameter()
n_groups   <- parameter()
max_degree <- parameter()
beta       <- parameter()
gamma      <- parameter(0.1)

dim(neighbors) <- c(max_degree, n)
neighbors <- parameter(type = "integer")

dim(mask)            <- c(max_degree, n)
dim(susceptibility)  <- n
dim(transmisibility) <- n
dim(S_ini)           <- n
dim(I_ini)           <- n
dim(G)               <- c(n_groups, n)
mask            <- parameter()
susceptibility  <- parameter()
transmisibility <- parameter()
S_ini           <- parameter()
I_ini           <- parameter()
G               <- parameter(type = "integer")

# Per-node state (internal — not exposed via the "state" filter we ask
# dust for, but present in the model so dynamics are per-agent). cumFOI
# accumulates the integrated force of infection experienced by each node
# (assuming perpetual susceptibility). We need it only so we can compute
# per-group mean and mean-of-squares aggregates (see cumFOI_sum,
# cumFOI_sumsq below) and get CV^2 of exposure heterogeneity without
# per-node output.
initial(S[]) <- S_ini[i]
initial(I[]) <- I_ini[i]
initial(R[]) <- 0
initial(C[]) <- 0
initial(cumFOI[]) <- 0
dim(S) <- n
dim(I) <- n
dim(R) <- n
dim(C) <- n
dim(cumFOI) <- n

# FOI via adjacency list — same as stoch_mod_adj.
contrib[, ] <- mask[i, j] * I[neighbors[i, j]] * transmisibility[neighbors[i, j]]
dim(contrib) <- c(max_degree, n)

foi[] <- beta / N_total * sum(contrib[, i])
dim(foi) <- n

N_total <- sum(S) + sum(I) + sum(R)

p_SI[] <- 1 - exp(-susceptibility[i] * foi[i] * dt)
p_IR   <- 1 - exp(-gamma  * dt)
dim(p_SI) <- n

n_SI[] <- Binomial(S[i], p_SI[i])
n_IR[] <- Binomial(I[i], p_IR)
dim(n_SI) <- n
dim(n_IR) <- n

update(S[]) <- S[i] - n_SI[i]
update(I[]) <- I[i] + n_SI[i] - n_IR[i]
update(R[]) <- R[i] + n_IR[i]
update(C[]) <- C[i] + n_SI[i]
update(cumFOI[]) <- cumFOI[i] + foi[i] * dt

# Group-weighted transitions (n_groups x n intermediates). Cost is O(n_groups
# * n) per step — small compared with FOI sums when n_groups is small.
GnSI[, ]  <- G[i, j] * n_SI[j]
GnIR[, ]  <- G[i, j] * n_IR[j]
G_Sini[, ] <- G[i, j] * S_ini[j]
G_Iini[, ] <- G[i, j] * I_ini[j]
dim(GnSI)   <- c(n_groups, n)
dim(GnIR)   <- c(n_groups, n)
dim(G_Sini) <- c(n_groups, n)
dim(G_Iini) <- c(n_groups, n)

dim(S_g) <- n_groups
dim(I_g) <- n_groups
dim(R_g) <- n_groups
dim(C_g) <- n_groups

initial(S_g[]) <- sum(G_Sini[i, ])
initial(I_g[]) <- sum(G_Iini[i, ])
initial(R_g[]) <- 0
initial(C_g[]) <- 0

update(S_g[]) <- S_g[i] - sum(GnSI[i, ])
update(I_g[]) <- I_g[i] + sum(GnSI[i, ]) - sum(GnIR[i, ])
update(R_g[]) <- R_g[i] + sum(GnIR[i, ])
update(C_g[]) <- C_g[i] + sum(GnSI[i, ])

# Per-group cumulative-FOI moments. We track sum_i G[g,i] * cumFOI[i]
# and sum_i G[g,i] * cumFOI[i]^2 as running state, incrementing from
# the current step's foi (which is computed before updates apply).
# For the sum, `new = old + dt * sum(G * foi)` is exact. For the sum
# of squares, `new = old + 2*dt*sum(G*cumFOI*foi) + dt^2*sum(G*foi^2)`
# — the algebra of (cumFOI + foi*dt)^2. Both use pre-step cumFOI (the
# current-step foi acts against the old cumFOI in that expansion).
Gfoi[, ]        <- G[i, j] * foi[j]
GcumFOI_foi[, ] <- G[i, j] * cumFOI[j] * foi[j]
Gfoi2[, ]       <- G[i, j] * foi[j] * foi[j]
dim(Gfoi)        <- c(n_groups, n)
dim(GcumFOI_foi) <- c(n_groups, n)
dim(Gfoi2)       <- c(n_groups, n)

dim(cumFOI_sum)   <- n_groups
dim(cumFOI_sumsq) <- n_groups
initial(cumFOI_sum[])   <- 0
initial(cumFOI_sumsq[]) <- 0
update(cumFOI_sum[])   <- cumFOI_sum[i]   + dt * sum(Gfoi[i, ])
update(cumFOI_sumsq[]) <- cumFOI_sumsq[i] +
                          2 * dt * sum(GcumFOI_foi[i, ]) +
                          dt * dt * sum(Gfoi2[i, ])
