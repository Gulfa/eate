# Determinsitic SIR model with n compartments
# Time varying beta
dim(beta_day) <- N_steps
beta_day[] <- user()
N_steps <- user()
interpolation_time[] <- user()
dim(interpolation_time) <- N_steps
beta <- interpolate(interpolation_time, beta_day, "linear")


import <- user(0)

# Main compartments updates
deriv(S[]) <- -n_SI[i] + n_RS_in[i] # n_import[i]
deriv(I[]) <-  n_SI[i] - n_IR[i]# + n_import[i]
deriv(R[]) <- n_IR[i]  - n_RS[i]
deriv(C[]) <- n_SI[i]
deriv(exposure[]) <- sum(lambda_ij[i,])


## Draws from binomial distributions for numbers changing between
## compartments:
n_SI[] <- S[i]*susceptibility[i]*sum(lambda_ij[i,])
n_IR[] <- I[i]*gamma
n_RS[] <- R[i]*waning
n_RS_in[] <- 0
n_RS_in[1] <- sum(n_RS[])

output(n_RS_in) <- TRUE
output(n_RS) <- TRUE
# Forece of infection
lambda_ij[,] <- beta * mixing_matrix[i,j]*I[j]/sum(N)*transmisibility[j]

#n_import[] <- rbinom(S[i], import/sum(S))

## Total population sizemet
N[] <- S[i] + I[i] + R[i]

## Initial states:
initial(S[]) <- S_ini[i]
initial(I[]) <- I_ini[i]
initial(R[]) <- 0
initial(C[]) <-  0
initial(exposure[]) <- 0
dim(n_SI) <-  n
dim(n_IR) <-n
dim(n_RS) <- n
dim(n_RS_in) <- n
dim(S) <- n
dim(I) <- n
dim(R) <- n
dim(N) <- n
dim(C) <- n

dim(exposure) <- n

dim(lambda_ij) <- c(n,n)
dim(beta_norm) <- n
dim(susceptibility) <- n
dim(transmisibility) <- n
output(n_SI) <- TRUE
# User defined parameters
gamma <- user(0.1)
n <- user(4)
S_ini[] <- user()
I_ini[] <- user()
mixing_matrix[,] <- user()
dim(mixing_matrix) <- c(n,n)
dim(S_ini) <- n
dim(I_ini) <- n
beta_norm[] <- user()
susceptibility[] <- user()
transmisibility[] <- user()
waning <- user(0)
              
