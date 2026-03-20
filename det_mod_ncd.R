

# Main compartments updates
deriv(S[]) <- -n_SR[i] # n_import[i]
deriv(R[]) <- n_SR[i]
deriv(C[]) <- n_SR[i]

hazzard_S[] <- n_SR[i]/S[i]
dim(hazzard_S) <- n
output(hazzard_S) <- TRUE

## Draws from binomial distributions for numbers changing between
## compartments:
n_SR[] <- S[i]*susceptibility[i]
output(n_SR) <- TRUE

## Total population size
N[] <- S[i] + R[i]

## Initial states:
initial(S[]) <- S_ini[i]
initial(R[]) <- 0
initial(C[]) <-  0

dim(n_SR) <- n
dim(S) <- n
dim(R) <- n
dim(N) <- n
dim(C) <- n

n <- user(4)
S_ini[] <- user()
susceptibility[] <- user()
dim(susceptibility) <- n
dim(S_ini) <- n

              
