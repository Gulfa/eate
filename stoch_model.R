library(adaptivetau)
library(data.table)

transitions_sir <- list(
  c(S1=-1, I1=1, C1=1), c(I1=-1, R1=1),
  c(S2=-1, I2=1, C2=1), c(I2=-1, R2=1)
)

transitions_linear <- list(
  c(S1=-1, C1=1),
  c(S2=-1, C2=1)
)

sir_rates <- function(x, p, t) {
  I_total <- x["I1"] + x["I2"]
  c(p$beta * x["S1"] * I_total / p$N,
    p$gamma * x["I1"],
    p$alpha * p$beta * x["S2"] * I_total / p$N,
    p$gamma * x["I2"])
}

linear_rates <- function(x, p, t) {
  c(p$beta * x["S1"],
    p$alpha * p$beta * x["S2"])
}

run_sir <- function(N_cont, N_vac, beta, alpha, gamma=1/7, I0_cont=1, I0_vac=0, t=100,
                    timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  params <- list(N=N_cont + N_vac, beta=beta, alpha=alpha, gamma=gamma)
  init <- c(S1=N_cont - I0_cont, I1=I0_cont, R1=0, C1=0,
            S2=N_vac  - I0_vac,  I2=I0_vac,  R2=0, C2=0)

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions_sir, sir_rates, params, tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

run_linear <- function(N_cont, N_vac, beta, alpha, t=100,
                       timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  params <- list(beta=beta, alpha=alpha)
  init <- c(S1=N_cont, C1=0, S2=N_vac, C2=0)

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions_linear, linear_rates, params, tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

run_stoch_cd <- function(mixing_matrix, beta, N, t, I_ini,
                         susceptibility=NULL, transmissibility=NULL,
                         gamma=1/3, timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n <- nrow(mixing_matrix)
  if (is.null(susceptibility))   susceptibility   <- rep(1, n)
  if (is.null(transmissibility)) transmissibility <- rep(1, n)
  N_total <- sum(N)

  S_names <- paste0("S", 1:n)
  I_names <- paste0("I", 1:n)
  R_names <- paste0("R", 1:n)
  C_names <- paste0("C", 1:n)

  transitions <- unlist(recursive=FALSE, lapply(1:n, function(i) {
    inf <- setNames(c(-1L,  1L, 1L), c(S_names[i], I_names[i], C_names[i]))
    rec <- setNames(c(-1L,  1L),     c(I_names[i], R_names[i]))
    list(inf, rec)
  }))

  rate_func <- local({
    mm  <- mixing_matrix
    sus <- susceptibility
    tr  <- transmissibility
    function(x, p, t) {
      I_vec <- x[I_names]
      foi   <- as.numeric(mm %*% (I_vec * tr)) / N_total
      rates <- numeric(2 * n)
      rates[seq(1, 2*n, 2)] <- pmax(0, sus * beta * x[S_names] * foi)
      rates[seq(2, 2*n, 2)] <- pmax(0, gamma * I_vec)
      rates
    }
  })

  init <- setNames(
    c(rbind(N - I_ini, I_ini, rep(0, n), rep(0, n))),
    c(rbind(S_names,   I_names, R_names,   C_names))
  )

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions, rate_func, list(), tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

run_stoch_linear <- function(beta, N, t, susceptibility=NULL,
                             timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  n <- length(N)
  if (is.null(susceptibility)) susceptibility <- rep(1, n)

  S_names <- paste0("S", 1:n)
  C_names <- paste0("C", 1:n)

  transitions <- lapply(1:n, function(i)
    setNames(c(-1L, 1L), c(S_names[i], C_names[i]))
  )

  rate_func <- local({
    sus <- susceptibility
    function(x, p, t) pmax(0, sus * beta * x[S_names])
  })

  init <- setNames(c(rbind(N, rep(0, n))), c(rbind(S_names, C_names)))

  res <- parallel::mclapply(1:n_sim, function(i) {
    out <- as.data.frame(ssa.adaptivetau(init, transitions, rate_func, list(), tf=t))
    regularise(out, timepoints) |> transform(sim=i)
  }, mc.cores=cores)
  rbindlist(res)
}

run_stoch_frailty_linear <- function(alpha, sd, beta=1, f=0.5, N=1000, t=100,
                                      n_frailty=100, vac_counts=NULL,
                                      timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  fr      <- get_frailty(sd=sd, n=n_frailty)
  frailty <- exp(2.5 * fr$x)
  n_total <- round(2 * N * fr$p)
  if (is.null(vac_counts)) vac_counts <- round(f * n_total)

  n_groups       <- 2 * n_frailty
  N_groups       <- c(n_total - vac_counts, vac_counts)
  susceptibility <- c(frailty, alpha * frailty)

  raw <- run_stoch_linear(beta=beta, N=N_groups, t=t, susceptibility=susceptibility,
                           timepoints=timepoints, n_sim=n_sim, cores=cores)

  N_vac   <- sum(vac_counts)
  N_unvac <- sum(n_total - vac_counts)
  raw[, vac   := rowSums(.SD), .SDcols=paste0("C", (n_frailty+1):n_groups)]
  raw[, unvac := rowSums(.SD), .SDcols=paste0("C", 1:n_frailty)]
  raw[, CRR   := (vac / N_vac) / (unvac / N_unvac)]
  raw[, .(time, sim, vac, unvac, CRR)]
}

run_stoch_frailty_cd <- function(alpha, sd, beta=1, R=NULL, f=0.5, N=1000, t=100,
                                  n_frailty=100, gamma=1/2, vac_counts=NULL,
                                  I_ini_total=1, timepoints=seq(0, t, 1), n_sim=100, cores=10) {
  if (!is.null(R)) beta <- get_beta(R, alpha, sd, f=f, N=N, n_frailty=n_frailty, gamma=gamma)

  fr      <- get_frailty(sd=sd, n=n_frailty)
  frailty <- exp(2.5 * fr$x)
  n_total <- round(2 * N * fr$p)
  if (is.null(vac_counts)) vac_counts <- round(f * n_total)

  n_groups       <- 2 * n_frailty
  N_groups       <- c(n_total - vac_counts, vac_counts)
  susceptibility <- c(frailty, alpha * frailty)
  mm             <- matrix(1, nrow=n_groups, ncol=n_groups) / n_groups

  I_ini <- rep(0L, n_groups)
  I_ini[which.max(N_groups[1:n_frailty])] <- I_ini_total

  raw <- run_stoch_cd(mm, beta=beta, N=N_groups, t=t, I_ini=I_ini,
                      susceptibility=susceptibility, gamma=gamma,
                      timepoints=timepoints, n_sim=n_sim, cores=cores)

  N_vac   <- sum(vac_counts)
  N_unvac <- sum(n_total - vac_counts)
  raw[, vac   := rowSums(.SD), .SDcols=paste0("C", (n_frailty+1):n_groups)]
  raw[, unvac := rowSums(.SD), .SDcols=paste0("C", 1:n_frailty)]
  raw[, CRR   := (vac / N_vac) / (unvac / N_unvac)]
  raw[, .(time, sim, vac, unvac, CRR)]
}

regularise <- function(df, timepoints) {
  df_reg <- dplyr::bind_rows(df, data.frame(time=timepoints)) |> dplyr::arrange(time)
  for (col in setdiff(colnames(df_reg), c("time", "sim"))) {
    df_reg[[col]] <- zoo::na.locf(df_reg[[col]], na.rm=FALSE)
  }
  dplyr::distinct(dplyr::filter(df_reg, time %in% timepoints))
}
