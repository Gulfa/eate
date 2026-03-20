




run_stacking <- function(init.values1, init.values2, transitions, params1, params2, rates1, rates2, T){

    n_species <- length(init.values1)
    species_names <- names(init.values1)
    n_reactions <- length(transitions)
    trans_matrix <- matrix(0, nrow = n_species, ncol = n_reactions)
    rownames(trans_matrix) <- species_names
    
    for (j in 1:length(transitions)) {
      trans <- transitions[[j]]
      if (is.null(names(trans))) {
        stop("Transitions in list format must be named vectors")
      }
      for (var_name in names(trans)) {
        var_idx <- which(species_names == var_name)
        if (length(var_idx) == 0) {
          stop(paste("Variable", var_name, "in transition not found in initial.values"))
        }
        trans_matrix[var_idx, j] <- trans[var_name]
      }
    }
    # transitions <- trans_matrix
    t <- 0
    ts <- c(t)
    

    Xs <- matrix(init.values1, nrow = 1)
    Zs <- matrix(init.values2, nrow = 1)
    X <- init.values1
    Z <- init.values2
 
    while (t < T) {
       
        rate_X <-rates1(X, params1, t)
        rate_Z <-rates1(Z, params2, t)
        
        lam_l <- pmax(rate_X, rate_Z)
        q <- c(0, cumsum(lam_l))
        lt <- sum(lam_l)
        if(lt<=0){
             ts <- c(ts, T)
            Xs <- rbind(Xs, X)
            Zs <- rbind(Zs, Z)
            break
        }
    e <- rexp(1, rate = 1)
    delta <- e/lt

    unif <- runif(1, 0, 1)
    mu <- 1e5

    for(i in 2:length(q)){
        if(q[i-1]/lt < unif & q[i]/lt > unif){
            mu <- i
            break
        }
    }
    
    if(unif < (q[mu -1] + rate_X[mu-1])/lt){
        X <- X + trans_matrix[,mu - 1]
    }
    if(unif < (q[mu -1] + rate_Z[mu-1])/lt){
        Z <- Z + trans_matrix[,mu -1]
    }
    t <- t + delta
    ts <- c(ts, t)
    Xs <- rbind(Xs, X)
    Zs <- rbind(Zs, Z)
    }
    return(list(X=as.data.frame(cbind(Xs, time=ts)),
                Z=as.data.frame(cbind(Zs, time=ts))))

}

test <- function(){

# Initial state (same for both)
initial.values <- c(N = 50)

# Transitions (rows are reactions, columns are state changes)
transitions <- list(c(N=+1), c(N=-1))
# Rate functions
rate_function_bd <- function(x, params, t) {
  c(
    params$beta * x["N"],    # Birth rate
    params$delta * x["N"]    # Death rate
  )
}

# Parameters
params1 <- list(beta = 0.03, delta = 0.02)  # Factual
params2 <- list(beta = 0.05, delta = 0.02)  # Counterfactual (higher birth rate)

# Run coupled simulation
result_bd <- run_stacking(  
    initial.values,
    initial.values,
    transitions,
    params1,
    params2,
    rate_function_bd,
    rate_function_bd,
    T = 10
)


N <- 300


res <- parallel::mclapply(1:N, function(i) {
  run_stacking(  
    initial.values,
    initial.values,
    transitions,
    params1,
    params2,
    rate_function_bd,
    rate_function_bd,
    T = 20
)
}, mc.cores = 8)
}




comb <- function(l){
    out <- data.frame()
    for(i in 1:length(l)){
        out <- rbind(out, l[[i]]$X %>% mutate(sim=i, type="X"))
        out <- rbind(out, l[[i]]$Z %>% mutate(sim=i, type="Z"))
    }
    return(out)
}
regularise <- function(df, timepoints){
    out <- data.frame()
    df_i_reg <- data.frame(time=timepoints)

    df_i_reg <- dplyr::bind_rows(df_i_reg, df) %>% arrange(time)
    
    for(col in colnames(df_i_reg)){
        if(!col %in% c("time", "sim")){
            df_i_reg[[col]] <- zoo::na.locf(df_i_reg[[col]], na.rm=FALSE)
        }

    }
    return(distinct(df_i_reg %>% filter(time %in% timepoints)))
}

compare_cf <- function(initial.values1,initial.values2, transitions, params1, params2, rate_function, T,  expr,N=10, n_cores=8, paramlist=FALSE){

    expr <- enquo(expr)
   
    if(paramlist){
        N_each <- N/length(params1)
        params1 <- rep(params1, N_each)
        params2 <- rep(params2, N_each)

    } else {
        
        params1 <- lapply(1:N, function(i) params1)#[[1]]
        params2 <- lapply(1:N, function(i) params2)#[[1]]
    
    }
    
      #  res <- parallel::mclapply(1:N, function(i) {
       #    # print(params1[[i]])
       # run_stacking(  
       #     initial.values1,
       #     initial.values2,
       #     transitions,
       ##     params2[[i]],
        #    rate_function,
        #    rate_function,
        #    T = T
        #)
        #}, mc.cores = n_cores)
        #res <- comb(res)
        #print("finished stacked sims")
        ind1 <- rbindlist(parallel::mclapply(1:N, function(i) regularise(as.data.frame(ssa.adaptivetau(initial.values1, transitions, rate_function, params1[[i]], tf=T)), seq(0, T, 0.1))%>%mutate(sim=i), mc.cores = n_cores))
        ind2 <- rbindlist(parallel::mclapply(1:N, function(i) regularise(as.data.frame(ssa.adaptivetau(initial.values2, transitions, rate_function, params2[[i]], tf=T)), seq(0, T, 0.1))%>%mutate(sim=i), mc.cores = n_cores))

    
    #print(regularise(res %>% filter(sim == i, type == "X"), seq(0,T,0.1)))
    # diff <- rbindlist(lapply(1:N, function(i) {
    #     res_x <- regularise(res %>% filter(sim == i, type == "X") %>%select(-sim), seq(0,T,0.1))
    #     res_z <- regularise(res %>% filter(sim == i, type == "Z") %>% select(-sim), seq(0,T,0.1))
    #     merged <- merge(res_x, res_z, by = "time", suffixes = c("_X", "_Z"))
        
    #     merged <- merged %>% mutate(res=!!expr)
    #     data.frame(time = merged$time, res = merged$res, sim = i) 
    # }))

    diff2 <- rbindlist(lapply(1:N, function(i) {
        ind1_i <- ind1 %>% filter(sim == i)
        ind2_i <- ind2 %>% filter(sim == i)
        merged <- merge(ind1_i, ind2_i, by = "time", suffixes = c("_X", "_Z")) %>% mutate(res=!!expr)
        data.frame(time = merged$time, res = merged$res, sim = i)
    }))

    tot <- rbind(#diff %>% mutate(type="stacked"),
                diff2 %>% mutate(type="ind"))
    return(list(diffs=tot, #stacked=res, 
                inds=rbind(ind1 %>% mutate(type="X"),
                     ind2 %>% mutate(type="Z"))))




}


def_run_test_linear <- function(n_cores=10){
    res <- compare_cf(
        c(N=100),
        c(N=100),
        list(c(N=1), c(N=-1)),
        list(birth=0.02, death=0.02),
        list(birth=0.03, death=0.03),
        function(x, p, t) {return(c(x*p$birth, x*p$death))},
        40,
        N_Z - N_X,
        N=200,
        n_cores=n_cores,
        paramlist=FALSE
    )
    test_sds(res)
}

test_sds <- function(res){
    
    ind1 <- res$inds %>% filter(type=="X")
    ind2 <- res$inds %>% filter(type=="Z")
    res <- res$stacked
    s1 <- sd(res %>% filter(type=="X") %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))# - res %>% filter(type=="Z") %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))
    s2 <- sd(ind1 %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))# - ind2 %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))
    print(glue::glue("Sds for process 1, {s1} and {s2}"))
    s3 <- sd(res %>% filter(type=="Z") %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))# - res %>% filter(type=="Z") %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))
    s4 <- d(ind2 %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))# - ind2 %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N))
    print(glue::glue("Sds for process 2, {s2} and {s3}"))
    s5 <- sd((res %>% filter(type=="X") %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N) - res %>% filter(type=="Z") %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N)))
    s6 <- sd((ind1 %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N) - ind2 %>% group_by(sim) %>% slice_tail()%>% ungroup()%>% pull(N)))
    print(glue::glue("Sds fordiff, {s5} and {s6}"))
}

test_nl <- function(){
        init.vals <-
            c(S1 = 1000,
                I1 = 3,
                R1 = 0,
                C1 = 0,
                S2 = 1000,
                I2 = 3,
                R2 = 0,
                C2 = 0
            ) # recovered (and immune) humans
        
        
        init.vals2 <-
            c(S1 = 999,
                I1 = 3,
                R1 = 0,
                C1 = 0,
                S2 = 1001,
                I2 = 3,
                R2 = 0,
                C2 = 0
            ) # recovered (and immune) humans
        
        
        
        transitions.vacSIR = list(
        c(S1 = -1, I1=1, C1=1),
        c(I1 = -1, R1=1),
        c(S2 = -1, I2=1, C2=1),
        c(I2 = -1, R2=1) 
        )

        params_nl = list(N=2006, beta=2, alpha=0.5)

    NonLinearRateF <- function(x, p, t) {
        return(c(p$beta*x["S1"]*(x["I1"] + x["I2"])/p$N, x["I1"], p$alpha*p$beta*x["S2"]*(x["I1"] + x["I2"])/p$N, x["I2"])
                ) # infection rate
    }


    sir_dif <- compare_cf(
        init.vals,
        init.vals2,
        transitions.vacSIR,
        
        list(params_nl, params_nl),
        list(params_nl, params_nl),
        NonLinearRateF,
        40,
        C2_Z / C1_X,
        N=50,
        n_cores=1,
        paramlist=TRUE
    )
}

# ggplot(tot, aes(x=time, y=diff, group=paste(type, sim), color=type) )+ geom_step() + theme_minimal()

# tot %>% group_by(type, time) %>% summarise(var=var(diff), m=mean(diff)) %>%
#     ggplot(aes(x=time, y=m, color=type)) + geom_line() + theme_minimal()

# ggplot(res %>% filter(sim %in% c(1:N)), aes(x=time, y=N, color=type, group=interaction(type, sim))) + geom_step() + theme_minimal()



# tot <- sir_dif$diffs

# ggplot(tot, aes(x=time, y=res, group=paste(type, sim), color=type) )+ geom_step() + theme_minimal()

# ggplot(tot %>% filter(time==40), aes(x=res, fill=type)) + geom_density(alpha=0.5) + theme_minimal() + xlim(0,5)

# tot %>% group_by(type, time) %>% summarise(var=var(diff), m=mean(diff)) %>%
#     ggplot(aes(x=time, y=m, color=type)) + geom_line() + theme_minimal()

# ggplot(res %>% filter(sim %in% c(1:N)), aes(x=time, y=N, color=type, group=interaction(type, sim))) + geom_step() + theme_minimal()




DE_sir <- function(N_cont, N_vac, I0_cont, I0_vac, T, param_list, n_sim=1, cores=10){
    
    state <- array(0, dim=c(N_cont + N_vac, 3))
    state[,1] <- 1
    suscpetibility <- c(rep(1, N_cont), rep(param_list$alpha, N_vac))
    state[1:I0_cont, 1] <- 0
    state[1:I0_cont, 2] <- 1

    with_vacs <- c()
    without_vacs <- c()



    res <- parallel::mclapply(1:n_sim, function(i){
        contacts <- generate_contacts(N_cont + N_vac, T=T, contact_rate=param_list$beta/(N_cont+N_vac))
        recovery_times <- rexp(I0_cont + I0_vac + length(contacts), rate=1)
        
        main_outcome <- run_sir_contacts(state, contacts, recovery_times, suscpetibility, param_list, T)

        with_vac <- 0
        without_vac <- 0
        for(j in 1:(N_cont + N_vac)){
            if(state[j,1]!=1) next
            new_suc <- copy(suscpetibility)
            if(j <= N_cont){
                new_suc[j] <- param_list$alpha
            } else {
                new_suc[j] <- 1
            }
            counterfactual_outcome <- run_sir_contacts(state, contacts, recovery_times, new_suc, param_list, T)
            main_j_outcome <- 1- main_outcome$last_state[j,1]
            counterfactual_j_outcome <- 1- counterfactual_outcome$last_state[j,1]


            if(j <= N_cont){
               with_vac <- with_vac + counterfactual_j_outcome
               without_vac <- without_vac + main_j_outcome
            } else {
                with_vac <- with_vac + main_j_outcome
                without_vac <- without_vac + counterfactual_j_outcome
            }
            
        }
        return(list(vaccinated=with_vac, unvaccinated=without_vac, tot_inf_main=N_cont + N_vac - I0_cont - I0_vac - sum(main_outcome$last_state[, 1]),
        tot_inf_unvac = N_cont - sum(main_outcome$last_state[I0_cont:N_cont, 1]),
        tot_inf_vac = N_vac - sum(main_outcome$last_state[(N_cont+1):(N_cont + N_vac), 1])) )}, mc.cores=cores)

    return(list(ratios = unlist(lapply(res, function(x) (x$vaccinated + 1e-9)/(x$unvaccinated + 1e-9))),

                diffs = unlist(lapply(res, function(x) x$vaccinated/(N_vac + N_cont - I0_vac - I0_cont) - x$unvaccinated/(N_vac + N_cont - I0_vac - I0_cont))),
                vacs = unlist(lapply(res, function(x) x$vaccinated)),
                unvacs = unlist(lapply(res, function(x) x$unvaccinated)),
                final_sizes = unlist(lapply(res, function(x) x$tot_inf_main)),
                final_sizes_unvac = unlist(lapply(res, function(x) x$tot_inf_unvac)),
                final_sizes_vac=unlist(lapply(res, function(x) x$tot_inf_vac))
                ))
                

}



DE_lin <- function(N_cont, N_vac, T, param_list, n_sim=1, cores=10){
    
    state <- rep(0, N_cont + N_vac)
    suscpetibility <- c(rep(1, N_cont), rep(param_list$alpha, N_vac))
    
    res <- parallel::mclapply(1:n_sim, function(i){
        events <- generate_linear_event_times(N_cont + N_vac, rate=param_list$beta, T=T)
                
        main_outcome <- run_events_linear(state, suscpetibility, events, T)
        with_vac <- 0
        without_vac <- 0
        for(j in 1:(N_cont + N_vac)){

            new_suc <- copy(suscpetibility)
            if(j <= N_cont){
                new_suc[j] <- param_list$alpha
            } else {
                new_suc[j] <- 1
            }
            counterfactual_outcome <- run_events_linear(state, new_suc, events, T)
            main_j_outcome <-  main_outcome$last_state[j]
            counterfactual_j_outcome <- counterfactual_outcome$last_state[j]
       
            if(j <= N_cont){
               with_vac <- with_vac + counterfactual_j_outcome
               without_vac <- without_vac + main_j_outcome
            } else {
                with_vac <- with_vac + main_j_outcome
                without_vac <- without_vac + counterfactual_j_outcome
            }
            
        }
        return(list(vaccinated=with_vac, unvaccinated=without_vac,tot_inf_unvac = sum(main_outcome$last_state[1:N_cont]), tot_inf_main=sum(main_outcome$last_state)))}, mc.cores=cores)

    return(list(ratios = unlist(lapply(res, function(x) (x$vaccinated + 1e-9)/(x$unvaccinated + 1e-9))),
                vacs = unlist(lapply(res, function(x) x$vaccinated)),
                unvacs= unlist(lapply(res, function(x) x$unvaccinated)),
                diffs = unlist(lapply(res, function(x) x$vaccinated/(N_vac + N_cont - I0_vac - I0_cont) - x$unvaccinated/(N_vac + N_cont - I0_vac - I0_cont))),
                final_sizes = unlist(lapply(res, function(x) x$tot_inf_main)),
                final_sizes_unvac = unlist(lapply(res, function(x) x$tot_inf_unvac))
                ))
                

}


DE_sir <- function(N_cont, N_vac, I0_cont, I0_vac, T, param_list, n_sim=1, cores=10){
    
    state <- array(0, dim=c(N_cont + N_vac, 3))
    state[,1] <- 1
    suscpetibility <- c(rep(1, N_cont), rep(param_list$alpha, N_vac))
    state[1:I0_cont, 1] <- 0
    state[1:I0_cont, 2] <- 1

    with_vacs <- c()
    without_vacs <- c()



    res <- parallel::mclapply(1:n_sim, function(i){
        contacts <- generate_contacts(N_cont + N_vac, T=T, contact_rate=param_list$beta/(N_cont+N_vac))
        recovery_times <- rexp(I0_cont + I0_vac + length(contacts), rate=1)
        
        main_outcome <- run_sir_contacts(state, contacts, recovery_times, suscpetibility, param_list, T)

        with_vac <- 0
        without_vac <- 0
        for(j in 1:(N_cont + N_vac)){
            if(state[j,1]!=1) next
            new_suc <- copy(suscpetibility)
            if(j <= N_cont){
                new_suc[j] <- param_list$alpha
            } else {
                new_suc[j] <- 1
            }
            counterfactual_outcome <- run_sir_contacts(state, contacts, recovery_times, new_suc, param_list, T)
            main_j_outcome <- 1- main_outcome$last_state[j,1]
            counterfactual_j_outcome <- 1- counterfactual_outcome$last_state[j,1]


            if(j <= N_cont){
               with_vac <- with_vac + counterfactual_j_outcome
               without_vac <- without_vac + main_j_outcome
            } else {
                with_vac <- with_vac + main_j_outcome
                without_vac <- without_vac + counterfactual_j_outcome
            }
            
        }
        return(list(vaccinated=with_vac, unvaccinated=without_vac, tot_inf_main=N_cont + N_vac - I0_cont - I0_vac - sum(main_outcome$last_state[, 1]),
        tot_inf_unvac = N_cont - sum(main_outcome$last_state[I0_cont:N_cont, 1]),
        tot_inf_vac = N_vac - sum(main_outcome$last_state[(N_cont+1):(N_cont + N_vac), 1])) )}, mc.cores=cores)

    return(list(ratios = unlist(lapply(res, function(x) (x$vaccinated + 1e-9)/(x$unvaccinated + 1e-9))),

                diffs = unlist(lapply(res, function(x) x$vaccinated/(N_vac + N_cont - I0_vac - I0_cont) - x$unvaccinated/(N_vac + N_cont - I0_vac - I0_cont))),
                vacs = unlist(lapply(res, function(x) x$vaccinated)),
                unvacs = unlist(lapply(res, function(x) x$unvaccinated)),
                final_sizes = unlist(lapply(res, function(x) x$tot_inf_main)),
                final_sizes_unvac = unlist(lapply(res, function(x) x$tot_inf_unvac)),
                final_sizes_vac=unlist(lapply(res, function(x) x$tot_inf_vac))
                ))
                

}



DE_lin <- function(N_cont, N_vac, T, param_list, n_sim=1, cores=10){
    
    state <- rep(0, N_cont + N_vac)
    suscpetibility <- c(rep(1, N_cont), rep(param_list$alpha, N_vac))
    
    res <- parallel::mclapply(1:n_sim, function(i){
        events <- generate_linear_event_times(N_cont + N_vac, rate=param_list$beta, T=T)
                
        main_outcome <- run_events_linear(state, suscpetibility, events, T)
        with_vac <- 0
        without_vac <- 0
        for(j in 1:(N_cont + N_vac)){

            new_suc <- copy(suscpetibility)
            if(j <= N_cont){
                new_suc[j] <- param_list$alpha
            } else {
                new_suc[j] <- 1
            }
            counterfactual_outcome <- run_events_linear(state, new_suc, events, T)
            main_j_outcome <-  main_outcome$last_state[j]
            counterfactual_j_outcome <- counterfactual_outcome$last_state[j]
       
            if(j <= N_cont){
               with_vac <- with_vac + counterfactual_j_outcome
               without_vac <- without_vac + main_j_outcome
            } else {
                with_vac <- with_vac + main_j_outcome
                without_vac <- without_vac + counterfactual_j_outcome
            }
            
        }
        return(list(vaccinated=with_vac, unvaccinated=without_vac,tot_inf_unvac = N_cont - sum(main_outcome$last_state[1:N_cont]), tot_inf_main=sum(main_outcome$last_state)))}, mc.cores=cores)

    return(list(ratios = unlist(lapply(res, function(x) (x$vaccinated + 1e-9)/(x$unvaccinated + 1e-9))),
                vacs = unlist(lapply(res, function(x) x$vaccinated)),
                unvacs= unlist(lapply(res, function(x) x$unvaccinated)),
                diffs = unlist(lapply(res, function(x) x$vaccinated/(N_vac + N_cont - I0_vac - I0_cont) - x$unvaccinated/(N_vac + N_cont - I0_vac - I0_cont))),
                final_sizes = unlist(lapply(res, function(x) x$tot_inf_main)),
                final_sizes_unvac = unlist(lapply(res, function(x) x$tot_inf_unvac))
                ))
                

}


TO_sir <- function(N_cont1, N_vac1,N_cont2, N_vac2, I0_cont, I0_vac, T, param_list, n_sim=1, cores=10){
    
    state <- array(0, dim=c(N_cont1 + N_vac1, 3))
    state[,1] <- 1
    suscpetibility1 <- c(rep(1, N_cont1), rep(param_list$alpha, N_vac1))
    suscpetibility2 <- c(rep(1, N_cont2), rep(param_list$alpha, N_vac2))
    state[1:I0_cont, 1] <- 0
    state[1:I0_cont, 2] <- 1

 
    res <- parallel::mclapply(1:n_sim, function(i){
        contacts <- generate_contacts(N_cont1 + N_vac1, T=T, contact_rate=param_list$beta/(N_cont+N_vac))
        recovery_times <- rexp(I0_cont + I0_vac + length(contacts), rate=1)
        
        outcome1 <- run_sir_contacts(state, contacts, recovery_times, suscpetibility1, param_list, T)
        outcome2 <- run_sir_contacts(state, contacts, recovery_times, suscpetibility2, param_list, T)


        tot_cases1 <- regularise(data.table::rbindlist(lapply(1:length(outcome1$states), function(i) data.frame(time=outcome1$time[[i]], tot1=sum(outcome1$states[[i]][,1]==0)))), seq(0.1, T, 0.1))
        tot_cases2 <- regularise(data.table::rbindlist(lapply(1:length(outcome2$states), function(i) data.frame(time=outcome2$time[[i]], tot2=sum(outcome2$states[[i]][,1]==0)))), seq(0.1, T, 0.1))


        return(cbind(tot_cases1, tot2=tot_cases2$tot2, diff=tot_cases1$tot1 - tot_cases2$tot2,ratio=tot_cases1$tot1/tot_cases2$tot2 ,sim=i))}, mc.cores=cores)

    return(data.table::rbindlist(res))
}



TO_lin <- function(N_cont1, N_vac1,N_cont2, N_vac2, I0_cont, I0_vac, T, param_list, n_sim=1, cores=10){
    
    state <- rep(0, N_cont1 + N_vac1)
    suscpetibility1 <- c(rep(1, N_cont1), rep(param_list$alpha, N_vac1))
    suscpetibility2 <- c(rep(1, N_cont2), rep(param_list$alpha, N_vac2))
    
        res <- parallel::mclapply(1:n_sim, function(i){
        events <- generate_linear_event_times(N_cont1 + N_vac1, rate=param_list$beta, T=T)
                
        outcome1 <- run_events_linear(state, suscpetibility1, events, T)
        outcome2<- run_events_linear(state, suscpetibility2, events, T)
        


        tot_cases1 <- regularise(data.table::rbindlist(lapply(1:length(outcome1$states), function(i) data.frame(time=outcome1$time[[i]], tot1=sum(outcome1$states[[i]])))), seq(0.1, T, 0.1))
        tot_cases2 <- regularise(data.table::rbindlist(lapply(1:length(outcome2$states), function(i) data.frame(time=outcome2$time[[i]], tot2=sum(outcome2$states[[i]])))), seq(0.1, T, 0.1))


        return(cbind(tot_cases1, tot2=tot_cases2$tot2, diff=tot_cases1$tot1 - tot_cases2$tot2,ratio=tot_cases1$tot1/tot_cases2$tot2, sim=i))}, mc.cores=cores)

    return(data.table::rbindlist(res))
}

TO_lin_indp <- function(N_cont1, N_vac1,N_cont2, N_vac2, I0_cont, I0_vac, T, param_list, n_sim=1, cores=10){
    
    state <- rep(0, N_cont1 + N_vac1)
    suscpetibility1 <- c(rep(1, N_cont1), rep(param_list$alpha, N_vac1))
    suscpetibility2 <- c(rep(1, N_cont2), rep(param_list$alpha, N_vac2))
    
    res <- parallel::mclapply(1:n_sim, function(i){
        events <- generate_linear_event_times(N_cont1 + N_vac1, rate=param_list$beta, T=T)
                
        outcome1 <- run_events_linear(state, suscpetibility1, events, T)
        events <- generate_linear_event_times(N_cont1 + N_vac1, rate=param_list$beta, T=T)
        outcome2<- run_events_linear(state, suscpetibility2, events, T)
        


        tot_cases1 <- regularise(data.table::rbindlist(lapply(1:length(outcome1$states), function(i) data.frame(time=outcome1$time[[i]], tot1=sum(outcome1$states[[i]])))), seq(0.1, T, 0.1))
        tot_cases2 <- regularise(data.table::rbindlist(lapply(1:length(outcome2$states), function(i) data.frame(time=outcome2$time[[i]], tot2=sum(outcome2$states[[i]])))), seq(0.1, T, 0.1))


        return(cbind(tot_cases1, tot2=tot_cases2$tot2, diff=tot_cases1$tot1 - tot_cases2$tot2,ratio=tot_cases1$tot2/tot_cases2$tot2, sim=i))}, mc.cores=cores)

    return(data.table::rbindlist(res))
}

TO_sir_indp <- function(N_cont1, N_vac1,N_cont2, N_vac2, I0_cont, I0_vac, T, param_list, n_sim=1, cores=10){
    
    state <- array(0, dim=c(N_cont1 + N_vac1, 3))
    state[,1] <- 1
    suscpetibility1 <- c(rep(1, N_cont1), rep(param_list$alpha, N_vac1))
    suscpetibility2 <- c(rep(1, N_cont2), rep(param_list$alpha, N_vac2))
    state[1:I0_cont, 1] <- 0
    state[1:I0_cont, 2] <- 1

    with_vacs <- c()
    without_vacs <- c()
    res <- parallel::mclapply(1:n_sim, function(i){
        contacts <- generate_contacts(N_cont1 + N_vac1, T=T, contact_rate=param_list$beta/(N_cont+N_vac))
        recovery_times <- rexp(I0_cont + I0_vac + length(contacts), rate=1)
        
        outcome1 <- run_sir_contacts(state, contacts, recovery_times, suscpetibility1, param_list, T)
        contacts <- generate_contacts(N_cont1 + N_vac1, T=T, contact_rate=param_list$beta/(N_cont+N_vac))
        recovery_times <- rexp(I0_cont + I0_vac + length(contacts), rate=1)
        
        outcome2 <- run_sir_contacts(state, contacts, recovery_times, suscpetibility2, param_list, T)


        tot_cases1 <- regularise(data.table::rbindlist(lapply(1:length(outcome1$states), function(i) data.frame(time=outcome1$time[[i]], tot1=sum(outcome1$states[[i]][,1]==0)))), seq(0.1, T, 0.1))
        tot_cases2 <- regularise(data.table::rbindlist(lapply(1:length(outcome2$states), function(i) data.frame(time=outcome2$time[[i]], tot2=sum(outcome2$states[[i]][,1]==0)))), seq(0.1, T, 0.1))


        return(cbind(tot_cases1, tot2=tot_cases2$tot2, diff=tot_cases1$tot1 - tot_cases2$tot2,ratio=tot_cases1$tot1/tot_cases2$tot2, sim=i))}, mc.cores=cores)

    return(data.table::rbindlist(res))
}


test_linear_implementation <- function(){
 
    cores <- 10
    n_sim <- 1000
    state <- rep(0,100)
    T <- 4
    suscpetibility <- c(rep(1, 50), rep(0.5, 50))
    
    res <- parallel::mclapply(1:n_sim, function(i){
        events <- generate_linear_event_times(100, rate=0.12, T=T)
        main_outcome <- run_events_linear(state, suscpetibility, events, T)
    }, mc.cores=cores)
    

    LinearRateF <- function(x, p, t) {
    return(c(p$beta*x["S1"], 0, p$alpha*p$beta*x["S2"], 0)
            ) # infection rate
    }
        init.vals <-
        c(S1 = 50,
            I1 = 0,
            R1 = 0,
            C1 = 0,
            N1 = 50,
            S2 = 50,
            I2 = 0,
            R2 = 0,
            C2 = 0,
            N2 = 50
          
        ) # recovered (and immune) humans
    transitions.vacSIR = list(
    c(S1 = -1, I1=1, C1=1),
    c(I1 = -1, R1=1),
    c(S2 = -1, I2=1, C2=1),
    c(I2 = -1, R2=1)
    )
    ind1 <- parallel::mclapply(1:n_sim, function(i) {
        a <- ssa.adaptivetau(init.vals, transitions.vacSIR,LinearRateF, list(beta=0.12, alpha=0.5), tf=T)
        tmp <- a[nrow(a), ]
        return(data.frame(sim=i,C1=tmp["C1"], C2=tmp["C2"]))
    }, mc.cores=cores)
    ind1 <- rbindlist(lapply(ind1, function(x) data.frame(x)))

    print(mean(unlist(lapply(res, function(x) sum(x$last_state[1:50])))))
    print(mean(ind1$C1))
    print(mean(unlist(lapply(res, function(x) sum(x$last_state[51:100])))))
    print(mean(ind1$C2))
}


test_sir_implementation <- function(){
 

    cores <- 10
    n_sim <- 2000
    T <- 4

    state <- array(0, dim=c(60, 3))
    state[,1] <- 1
    suscpetibility <- c(rep(1, 30), rep(0.5, 30))
    state[1:2, 1] <- 0
    state[1:2, 2] <- 1



    res <- parallel::mclapply(1:n_sim, function(i){
        contacts <- generate_contacts(60, T=T, contact_rate=2/60)
        recovery_times <- rexp(2 + length(contacts), rate=1)

        return(run_sir_contacts(state, contacts, recovery_times, suscpetibility, param_list, T))
    }, mc.cores=cores)


    NonLinearRateF <- function(x, p, t) {
    return(c(p$beta*x["S1"]*(x["I1"] + x["I2"])/p$N, x["I1"], p$alpha*p$beta*x["S2"]*(x["I1"] + x["I2"])/p$N, x["I2"])
            ) # infection rate
    }
        init.vals <-
        c(S1 = 28,
            I1 = 2,
            R1 = 0,
            C1 = 0,
            N1 = 30,
            S2 = 30,
            I2 = 0,
            R2 = 0,
            C2 = 0,
            N2 = 30
          
        ) # recovered (and immune) humans
    transitions.vacSIR = list(
    c(S1 = -1, I1=1, C1=1),
    c(I1 = -1, R1=1),
    c(S2 = -1, I2=1, C2=1),
    c(I2 = -1, R2=1)
    )
    ind1 <- parallel::mclapply(1:n_sim, function(i) {
        a <- ssa.adaptivetau(init.vals, transitions.vacSIR,NonLinearRateF, list(beta=2,N=60, alpha=0.5), tf=T)
        tmp <- a[nrow(a), ]
        return(data.frame(sim=i,C1=tmp["C1"], C2=tmp["C2"]))
    }, mc.cores=cores)
    ind1 <- rbindlist(lapply(ind1, function(x) data.frame(x)))

    print(mean(unlist(lapply(res, function(x) 28 - sum(x$last_state[3:30,1] )))))
    print(mean(ind1$C1))
    print(mean(unlist(lapply(res, function(x) 30 - sum(x$last_state[31:60,1] )))))
    print(mean(ind1$C2))
}

run_sir_contacts <- function(init_state, contacts, recovery_times, suscpetibility, param_list, T){

    state <- init_state
    rec_time_counter <- 1
    rec_events <- list()
    for(i in which(state[,2]==1)){
        rec_events <- append(rec_events, list(list(type="recovery", which=i, t=recovery_times[rec_time_counter])))
        rec_time_counter <- rec_time_counter + 1
    }
    rec_time_counter <- rec_time_counter - 1
    events <-append(contacts, rec_events)
    sorted_events <- events[order(sapply(events, function(x) x$t))]

    states <- list(state)
    t <- 0
    times <- c(0)


   while(length(sorted_events) > 0 & t < T){
        event <- sorted_events[[1]]
        
        change <- FALSE
        if(event$type=="recovery"){
        
            
            new_state <- copy(state)
            new_state[event$which, 2] <- 0
            new_state[event$which, 3] <- 1
            change <-TRUE
        }
        if(event$type=="contact"){
            rec_time_counter <- rec_time_counter + 1
            if(state[event$from, 2]!= 1 | state[event$to, 1] != 1){
                sorted_events <- sorted_events[-1]
                next
            }
         
            new_state <- copy(state)

            if(event$inf_rand > suscpetibility[event$to]){
                sorted_events <- sorted_events[-1]
                next
            }

            new_state[event$to, 2] <- 1
            new_state[event$to, 1] <- 0
            sorted_events <- append(sorted_events, list(list(type="recovery", which=event$to, t=event$t + recovery_times[rec_time_counter] )))
            sorted_events <- sorted_events[order(sapply(sorted_events, function(x) x$t))]
            change <- TRUE
        }

        if(change){
            t <- event$t
            
            if(t < T){
                states <- append(states, list(new_state))
                state <- new_state
                times <- c(times, t)
                
            }
        }
        sorted_events <- sorted_events[-1]
        
    }
     states <- append(states, list(new_state))
     times <- c(times, T)
     
    return(list(states=states, time=times, last_state=state) )
}



generate_linear_event_times <- function(N, rate, T){
  
    events <- list()
    for(i in 1:N){
          t <- 0
    while(t < T){
        new_t <- rexp(1, rate=rate)
        if(t + new_t < T){
            events <- append(events, list(list(t=t+new_t, who=i, unif=runif(1)) ) )
        }
        t <- t + new_t
    }
    }
    return(events)
}

run_events_linear <- function(state, suscpetibility, events, T){
    sorted_events <- events[order(sapply(events, function(x) x$t))]
    states <- list(state)
    t <- 0
    times <- c(0)
    new_state <- copy(state)

    while(length(sorted_events) > 0 & t < T){
        event <- sorted_events[[1]]

        new_state <- copy(state)
        if(state[event$who] == 1){
            sorted_events <- sorted_events[-1]
            next
        }
        if(event$unif < suscpetibility[event$who]){
            new_state[event$who] <- 1
        }
        t <- event$t
        if(t < T){
            states <- append(states, list(new_state))
            state <- new_state
            times <- c(times, t)
        }
        sorted_events <- sorted_events[-1]

    }
     states <- append(states, list(new_state))
     times <- c(times, T)
    return(list(states=states, time=times, last_state=state) )
}

generate_contacts <- function(N,T, contact_rate){
    contacts <- list()
    for(i in 1:N){
        for(j in 1:N){
            
            if(i == j) next
            t <- 0
            while(t < T){
                new_t <- rexp(1, rate=contact_rate)
                if(t + new_t < T){
                    contacts <- append(contacts, list(list(type="contact",from=i,to=j, t=t + new_t, inf_rand=runif(1)) ))
                }
                t <- t + new_t
            }
        }
    }
    
    return(contacts)
}
