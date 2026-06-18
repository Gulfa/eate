

esteimate_VE_from_fs <- function(generator,eate_func, t=8, N_vac=100, X_vac=50, N_cont=100, X_cont=50, n_cores=10, beta_ini=0.1, alpha_ini=0.1, n=5000, burn_in=200, scale=c(0.05, 0.05)) {
    params <- fit_mod(generator, beta_ini=beta_ini, alpha_ini=alpha_ini, X_vac=X_vac, X_cont=X_cont, n=n, burn_in=burn_in, scale=scale)

    eate_func <- purrr::partial(eate_func, vac_frac=N_vac/(N_vac+N_cont), t=t, N=N_vac+N_cont)
    VE <- calc_ve_from_pars(eate_func, params, n_cores=n_cores)


    return(list(params=params, VE=VE))
}


# ---------------------------------------------------------------------------
# Config-driven model builders. Each takes a plain-data config and returns
# list(generator, eate_func). Shared state (c_ij for network, frailty bins
# for frailty) is materialised here, so configs stay small and serialisable
# for slurm. The EATE function always samples its own vaccine allocation
# internally — only the contact structure is shared with the fit generator.
# ---------------------------------------------------------------------------

build_network_mod <- function(cfg) {
    N_total <- sum(cfg$N)
    n_vac   <- cfg$N[2]

    set.seed(cfg$network_seed)
    c_ij <- get_conact_matrix_pl(N_total, alpha=cfg$pl_alpha, mean_k=cfg$mean_k)
    set.seed(NULL)

    set.seed(cfg$allocation_seed)
    vac <- sample(seq_len(N_total), n_vac)
    set.seed(NULL)

    generator <- purrr::partial(run_stoch_network,
                                c_ij=c_ij, vac=vac, N=N_total, k_mean=cfg$mean_k,
                                gamma=cfg$gamma, t=cfg$t, dt=cfg$dt, I_ini=cfg$I_ini,
                                method="dust", timepoints=cfg$timepoints,
                                n_sim=cfg$n_sim, cores=cfg$cores)

    # get_stoch_eate_network uses `f` instead of `vac_frac`; adapt the call
    # so the outer purrr::partial in esteimate_VE_from_fs (which binds
    # vac_frac/t/N) still works. n_rep stochastic replicates per allocation
    # replace the N+1 deterministic perturbation runs.
    eate_func <- function(beta, susceptibility, vac_frac, t, N, ...) {
        get_stoch_eate_network(beta=beta, susceptibility=susceptibility,
                               f=vac_frac, t=t, N=N, c_ij=c_ij,
                               k_mean=cfg$mean_k, init_I=cfg$I_ini,
                               gamma=cfg$gamma,
                               n_vac=cfg$eate_n_vac, n_rep=cfg$eate_n_rep,
                               mc.cores=cfg$cores)
    }
    list(generator=generator, eate_func=eate_func)
}

build_frailty_mod <- function(cfg) {
    # run_stoch_frailty_cd / get_frailty_eate use a "N-per-group" convention:
    # internally they build a population of size 2*N (vac + unvac). To match
    # the data total sum(cfg$N), we pass N = sum(cfg$N)/2.
    N_per_group <- sum(cfg$N) / 2

    sd_pop  <- if (cfg$sd > 0) cfg$sd else cfg$sd_trans
    fr      <- get_frailty(sd=sd_pop, n=cfg$n_frailty)
    n_total <- round(2 * N_per_group * fr$p)   # same formula the simulator uses
    bin     <- rep(seq_along(n_total), n_total)

    # cfg$N[2] is the number vaccinated in the data; tabulate that many
    # bin-labels into per-bin vaccinated counts.
    set.seed(cfg$allocation_seed)
    vacc <- tabulate(bin[sample(length(bin), cfg$N[2])], nbins=length(n_total))
    set.seed(NULL)

    generator <- purrr::partial(run_stoch_frailty_cd,
                                sd=cfg$sd, sd_trans=cfg$sd_trans,
                                I_ini_total=sum(cfg$I_ini), N=N_per_group,
                                gamma=cfg$gamma,
                                t=cfg$t, dt=cfg$dt, vac_counts=vacc,
                                n_frailty=cfg$n_frailty, timepoints=cfg$timepoints,
                                n_sim=cfg$n_sim, cores=cfg$cores, method="dust")

    # esteimate_VE_from_fs binds N=N_vac+N_cont (data total) to the eate_func.
    # The frailty EATE function uses the same 2*N convention as the simulator,
    # so override to N/2 inside. n_rep stochastic replicates replace the
    # deterministic frozen-field integration; I_ini_total is shared with the
    # generator so initial-seed dynamics match.
    eate_func <- function(beta, susceptibility, vac_frac, t, N, ...) {
        get_stoch_eate_frailty(alpha=susceptibility[2],
                               sd=cfg$sd, sd_trans=cfg$sd_trans,
                               beta=beta, f=vac_frac, N=N/2, t=t,
                               n_frailty=cfg$n_frailty,
                               gamma=cfg$gamma,
                               I_ini_total=sum(cfg$I_ini),
                               n_vac=cfg$eate_n_vac, n_rep=cfg$eate_n_rep,
                               mc.cores=cfg$cores)
    }
    list(generator=generator, eate_func=eate_func)
}


run_from_config <- function(config){
    model_type <- if (is.null(config$model_type)) "legacy" else config$model_type
    mod <- switch(model_type,
                  network = build_network_mod(config),
                  frailty = build_frailty_mod(config),
                  legacy  = list(generator=config$generator, eate_func=config$eate_func))

    res <- esteimate_VE_from_fs(generator=mod$generator, eate_func=mod$eate_func, t=config$t, N_vac=config$N_vac, X_vac=config$X_vac, N_cont=config$N_cont, X_cont=config$X_cont, n_cores=config$n_cores,
             beta_ini=config$beta_ini, alpha_ini=config$alpha_ini, n=config$n, burn_in=config$burn_in, scale=config$scale)
    res$name <- config$name
   return(res)

}


get_frailty_mod <- function(N=c(10,10), t=8, dt=0.1, timepoints=c(1,2,3), n_sim=10, cores=10, I_ini=c(10,10), allocation_seed=1, sd=0, sd_trans=0.3, n_frailty=10) {
    # WIP — body commented out to keep the file parseable.
    # for()
     sd_pop <- if (sd > 0) sd else sd_trans
     fr     <- get_frailty(sd=sd_pop, n=n_frailty)
     n_total <- round(sum(N) * fr$p)
    bin  <- rep(seq_along(n_total), n_total)
    set.seed(allocation_seed)
    vacc <- tabulate(bin[sample(length(bin), N[1])], nbins = length(n_total))
    set.seed(NULL)

    print(vacc)
    return(purrr::partial(run_stoch_frailty_cd, sd, sd_trans=sd_trans, I_ini_total=sum(I_ini), N=sum(N), t=t, dt=dt,vac_counts=vacc,n_frailty=n_frailty,
                             timepoints=timepoints, n_sim=n_sim, cores=cores, method="dust"))

}


get_network_mod <- function(N=c(100,100), pl_alpha=3, mean_k=6, t=8, dt=0.1,
                            timepoints=c(1,2,3), n_sim=10, cores=10, I_ini=2,
                            gamma=1/3, network_seed=1, allocation_seed=1) {
    N_total <- sum(N)
    n_vac   <- N[2]

    set.seed(network_seed)
    c_ij <- get_conact_matrix_pl(N_total, alpha=pl_alpha, mean_k=mean_k)
    set.seed(NULL)

    set.seed(allocation_seed)
    vac <- sample(seq_len(N_total), n_vac)
    set.seed(NULL)

    return(purrr::partial(run_stoch_network, c_ij=c_ij, vac=vac, N=N_total,
                          k_mean=mean_k, gamma=gamma, t=t, dt=dt, I_ini=I_ini,method="dust",    
                          timepoints=timepoints, n_sim=n_sim, cores=cores))
}


# Interactive scratch — guarded so the file can be sourced from Rscript
# (run_array.R) without trying to execute these calls.
if (FALSE) {

N_cont <- 100
X_cont <- 50
N_vac <- 100
X_vac <- 25
t <- 8



r <- run_from_config(list(name="linear", generator= purrr::partial(run_stoch_linear_dust, N=c(N_cont, N_vac), t=t, dt=0.01,
                             timepoints=seq(1, t, 1), n_sim=1000, cores=n_cores),
                              eate_func=get_EATE_linear, t=t, N_vac=N_vac,
                                X_vac=X_vac, N_cont=N_cont, X_cont=X_cont, n_cores=n_cores, beta_ini=0.05, alpha_ini=0.5, n=500, burn_in=500, scale=c(0.02, 0.05)))



r <- run_from_config(list(name="sir", generator= purrr::partial(run_stoch_cd_dust, matrix(rep(1, 4), nrow=2), I_ini=c(10,10), N=c(N_cont, N_vac), t=t, dt=0.01,
                             timepoints=seq(1, t, 1), n_sim=1000, cores=n_cores),
                              eate_func=get_EATE_linear, t=t, N_vac=N_vac,
                            X_vac=X_vac, N_cont=N_cont, X_cont=X_cont, n_cores=n_cores, beta_ini=1.5, alpha_ini=0.5, n=50, burn_in=50, scale=c(0.1, 0.05)))
r <- run_from_config(list(name="sir_frailty", generator= get_frailty_mod(N=c(N_cont, N_vac), t=t, dt=0.01, timepoints=seq(1, t, 1), n_sim=1000, cores=n_cores, I_ini=c(10,10), allocation_seed=1, sd=0, sd_trans=0.3, n_frailty=10),
    eate_func=get_EATE_linear, t=t, N_vac=N_vac,
     X_vac=X_vac, N_cont=N_cont, X_cont=X_cont, n_cores=n_cores, beta_ini=1.5, alpha_ini=0.5, n=50, burn_in=50, scale=c(0.1, 0.05)))

r <- run_from_config(list(name="sir_network", generator= get_network_mod(N=c(N_cont, N_vac), pl_alpha=3, mean_k=6, t=t, dt=0.01, timepoints=seq(1, t, 1), n_sim=1000, cores=n_cores, I_ini=2, network_seed=1, allocation_seed=1),
    eate_func=get_eate_network, t=t, N_vac=N_vac,
     X_vac=X_vac, N_cont=N_cont, X_cont=X_cont, n_cores=n_cores, beta_ini=1.5, alpha_ini=0.5, n=50, burn_in=50, scale=c(0.1, 0.05)))


# ---------------------------------------------------------------------------
# Pure-data config style for slurm: each list element is fully self-describing,
# c_ij / frailty bins are materialised inside run_from_config on the worker.
# ---------------------------------------------------------------------------

configs <- list(
    list(name="sir_frailty", model_type="frailty",
         N=c(N_cont, N_vac), t=t, dt=0.01,
         timepoints=seq(1, t, 1), n_sim=1000, cores=10,
         I_ini=c(10,10), allocation_seed=1,
         sd=0, sd_trans=0.3, n_frailty=10, gamma=1/2,
         eate_method="full", eate_slowdown=1, eate_n_vac=10,
         N_vac=N_vac, X_vac=X_vac, N_cont=N_cont, X_cont=X_cont,
         n_cores=10, beta_ini=1.5, alpha_ini=0.5, n=50, burn_in=50, scale=c(0.1, 0.05)),

    list(name="sir_network", model_type="network",
         N=c(N_cont, N_vac), pl_alpha=3, mean_k=6,
         t=t, dt=0.01, timepoints=seq(1, t, 1), n_sim=1000, cores=10,
         I_ini=2, gamma=1/3, network_seed=1, allocation_seed=1,
         eate_method="full", eate_slowdown=1, eate_n_vac=10,
         N_vac=N_vac, X_vac=X_vac, N_cont=N_cont, X_cont=X_cont,
         n_cores=10, beta_ini=1.5, alpha_ini=0.5, n=50, burn_in=50, scale=c(0.1, 0.05))
)

# Single-node parallel:
# results <- parallel::mclapply(configs, run_from_config, mc.cores=length(configs))
# Slurm: one job per config — write configs to RDS, dispatch one task per index.

}  # end scratch guard

