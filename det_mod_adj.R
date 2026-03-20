# Sparse adjacency-list SIR model for odin2.
#
# odin2 supports  result[] <- sum(matrix[, i])  where i is the LHS loop
# index — this was not possible in odin 1.x, which is why this model
# requires odin2.
#
# Network representation (precomputed by contact_matrix_to_adj() in model.R):
#   neighbors[max_degree, n]  — integer neighbor indices, padded with 1
#   mask[max_degree, n]       — 1.0 for real neighbors, 0.0 for padding
#
# FOI per node is O(max_degree) instead of O(n), giving O(edges) total.

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

# Step 1: per-slot infectivity contribution (max_degree x n)
contrib[,] <- mask[i, j] * I[neighbors[i, j]] * transmisibility[neighbors[i, j]]
dim(contrib) <- c(max_degree, n)

# Step 2: sum over neighbor slots — valid in odin2
foi[] <- beta / N_total * sum(contrib[, i])
dim(foi) <- n

N_total <- sum(S) + sum(I) + sum(R)

deriv(S[]) <- -S[i] * susceptibility[i] * foi[i] - R[i] * waning
deriv(I[]) <-  S[i] * susceptibility[i] * foi[i] - I[i] * gamma
deriv(R[]) <-  I[i] * gamma - R[i] * waning
deriv(C[]) <-  S[i] * susceptibility[i] * foi[i]
