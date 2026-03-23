# Try updating packages
#install.packages(c("odin", "glue"))
# If that doesn't work, try installing from GitHub
#remotes::install_github("mrc-ide/odin")

#odin_model <- odin.dust::odin_dust("odinmodel.R")

det_model_cd  <- odin::odin("det_mod_cd.R")
det_model_ncd <- odin::odin("det_mod_ncd.R")
det_model_adj <- odin2::odin("det_mod_adj.R")

#det_model_hazard <- odin::odin("det_hazard.R")
#det_non_linear <- odin::odin("simple_non_linear.R")


run_hazard_model <- function(hazard1, hazard2,t, N=1e4){

  model <- det_model_hazard$new(user=list(hazard1_day=hazard1, hazard2_day=hazard2, N_steps=t, S_ini=rep(N, 2), interpolation_time=1:t))
  res <- model$run(1:t)
  C1 <- N - res[,2]
  C2 <- N - res[,3]
  inc1 <- c(0, C1[2:length(C1)] - C1[1:(length(C1) -1)])
  inc2 <- c(0, C2[2:length(C2)] - C2[1:(length(C2) -1)])
   d <- as.data.frame(res)
  d$inc1 <- inc1
  d$inc2 <- inc2
  d$RR <- d$inc2/d$inc1
  d$cum_I_RR <- cumsum(d$inc2)/cumsum(d$inc1)

  return(d)
}

library(odin)

simple_model <- odin::odin({
  deriv(y) <- -0.5 * y
  initial(y) <- 1
})

mod <- simple_model()

run_frailty <- function(lambda, alpha, sd, N=1000, t=5000, n_frailty=200){
 

    if(sd>0){
           f <- get_frailty(sd=sd, n=n_frailty)
    frailty <- exp(log(lambda) + 2.5*f$x)
    
    params <- list(
        n=2*n_frailty,
        S_ini=N*rep(f$p,2),#*(f$x[2] - f$x[1]),
        susceptibility=c(frailty, alpha*frailty)
    )}
    else{
        params <- list(

        n=2,
        S_ini=c(N, N),
        susceptibility=c(lambda, lambda*alpha))
    }
    model <-det_model_ncd$new(user=params)
    res <- model$run(0:t)
    exp <- res %>% as.data.frame() %>%select(starts_with("C"))
    if(sd>0){
        df <- data.frame(t=res[,1], vac=rowSums(exp[, (params$n/2+1):(params$n)]), unvac=rowSums(exp[, 1:(params$n/2)]))
        
    }else{
        df <- data.frame(t=res[,1], vac=exp[,2], unvac=exp[, 1])
    }
    #browser()
  
      df$inc_unvac <- rowSums(res[, paste("n_SR[", 1:(params$n/2), "]", sep="")])
      df$inc_vac <- rowSums(res[, paste("n_SR[", (params$n/2+1):(params$n), "]", sep="")]) 
     df$HRR <- (df$inc_vac / (N - df$vac))/(df$inc_unvac / (N - df$unvac))
    #last_HRR_i <- max(which((N-df$vac)> 1e-4 & (N-df$unvac)> 1e-4 & df$inc_unvac>1e-7 & df$inc_vac>1e-7))
    #df$HRR[last_HRR_i:nrow(df)] <- df$HRR[last_HRR_i-1]

    strats <- list()
    for(i in 1:n_frailty){
        df_strat <- data.frame(t=res[,1], frailty_group=as.character(i))
        df_strat[,"vac"] <- res[, paste("C[", n_frailty + i, "]", sep="")] / params$S_ini[n_frailty + i]
        df_strat[, "unvac"] <- res[, paste("C[", i, "]", sep="")] / params$S_ini[i]
        df_strat$haz_unvac <- res[, paste("hazzard_S[", i, "]", sep="")]
        df_strat$haz_vac <- res[, paste("hazzard_S[", n_frailty + i, "]", sep="")]

        df_strat$HRR <- (df_strat$haz_vac)/(df_strat$haz_unvac)
        #last_HRR_i <- max(which((df_strat$vac< 0.99 & df_strat$unvac < 0.99 & (df_strat$inc_unvac)>1e-7 & df_strat$inc_vac>1e-7)))
        #print(last_HRR_i)
        df_strat$HRR[last_HRR_i:nrow(df_strat)] <- df_strat$HRR[last_HRR_i-1]
    
        df_strat$CRR  <- (df_strat$vac/df_strat$unvac)
         
      strats[[i]] <- df_strat



    }
    

    return(list(sum=df %>% mutate(CRR=vac/unvac, VE=1-CRR), strats=rbindlist(strats)))
}


get_conact_matrix_pl <- function(N, alpha, mean_k=6){


  
  propensities <- Pareto::rPareto(N, alpha=alpha, t=1)
  propensities <- propensities / mean(propensities)
  contacts <- matrix(0, nrow=N, ncol=N)
  for(i in 1:N){
      for(j in 1:N){
        
          contacts[i,j] <- rbinom(1, 1, min(mean_k*propensities[i]*propensities[j]/N, 1))
      }
  } 

  return(contacts)
}


run_mean_field <- function(beta=1, N=100, pl_alpha=3, alpha=1, t=100,vac_frac=0.5,vac=NULL, gamma=1/3, c_ij=NULL, k_mean=6){
  
  if(is.null(c_ij)){
      c_ij <- get_conact_matrix_pl(N, pl_alpha, mean_k=k_mean)
  }
  if(is.null(vac) ){
      vac <- sample(1:N, vac_frac*N)
  }
  susept <- rep(1, N)
  susept[vac] <- alpha
  non_vac <- setdiff(1:N, vac)
  res <- run_det_cd(c_ij, rep(N*beta/k_mean, t), rep(1,N), t, I_ini=c(1,1,rep(0, N-2)), gamma=gamma, beta_norm=rep(1, N), susceptibility=susept, sparse=TRUE)

  sum <- data.frame(t=res$full_results[["t"]], vac=rowSums(res$full_results[, paste("C[", vac, "]", sep="")]), unvac=rowSums(res$full_results[, paste("C[", non_vac, "]", sep="")]))
  sum$CRR <- (sum$vac/(vac_frac*N))/(sum$unvac/((1-vac_frac)*N))
sum  
  return(list(sum=sum, full=res$full_results))
}


get_eate_frailty <- function(alpha, sd,beta=1, R=NULL, f=0.5, N=1000, t=100, n_frailty=100){
    
  full_res <- run_frailty_cd(alpha, sd, N=N,beta=beta, R=R, t=t,f=f, n_frailty=n_frailty)
  full_eff <- 0
  for(k in 1:n_frailty){
      res_k <- run_frailty_cd(alpha, sd, N=N,beta=beta, R=R, t=t,f=f, n_frailty=n_frailty, k=k)
      res_mk <- run_frailty_cd(alpha, sd, N=N,beta=beta, R=R, t=t, f=f, n_frailty=n_frailty, k=-k)
      N_unvac  <- full_res$full[, paste("S[", k, "]", sep="")][1] 
      N_vac  <- full_res$full[, paste("S[", n_frailty + k, "]", sep="")][1]
      if(all(res_k==0) | all(res_mk==0)){
          next
      }
      eff <- res_k / (full_res$full[ paste("C[", k, "]", sep="")]/N_unvac)
      eff2 <-  (full_res$full[ paste("C[", n_frailty + k, "]", sep="")]/N_vac)/(res_mk)
      full_eff <- full_eff + eff*N_unvac + eff2 * N_vac
  }
  full_eff <- as.numeric(unlist(full_eff / (2*N)))

  return(full_res$sum %>% mutate(EATE=full_eff))

}

get_eate_network <- function(alpha=0.5,beta=1, R=NULL, f=0.5, N=200, t=15, pl_alpha=3, c_ij=NULL, n_vac=10){
    
  if(is.null(c_ij)){
      c_ij <- get_conact_matrix_pl(N, pl_alpha)
  }
  
  res <- list()

  run_vac <- function(){
    
     vac <- sample(1:N, f*N)
    full_res <- run_mean_field(beta=beta, N=N, alpha=alpha, t=t, vac_frac=f, gamma=1, c_ij=c_ij, vac=vac)
    denom <- rep(0, t)
    num <- rep(0, t)
    for(k in 1:N){

      if(k %in% vac){
          num <- num + full_res$full[, paste("C[", k, "]", sep="")]
          vac_to_run <- vac[vac!=k]
          res_mk <- run_mean_field(beta=beta, N=N, alpha=alpha, t=t, vac_frac=f, gamma=1, c_ij=c_ij, vac=vac_to_run)
          denom <- denom + res_mk$full[, paste("C[", k, "]", sep="")]
      }else{
          denom <- denom + full_res$full[, paste("C[", k, "]", sep="")]
          vac_to_run <- c(vac, k)
          res_k <- run_mean_field(beta=beta, N=N, alpha=alpha, t=t, vac_frac=f, gamma=1, c_ij=c_ij, vac=vac_to_run)
          num <- num + res_k$full[, paste("C[", k, "]", sep="")]
        }
      }
    return(data.frame(t=1:t, eate=num/denom, CRR=full_res$sum$CRR, sim=runif(1)))
  
  }

  res <- parallel::mclapply(1:n_vac, function(i) run_vac(), mc.cores=10)

  return(rbindlist(res))

}


run_frailty_cd <- function(alpha, sd,beta=1, R=NULL, f=0.5, N=1000, t=100, n_frailty=100, k=NULL, gamma=1/2){
 

    if(sd>0){
          fr <- get_frailty(sd=sd, n=n_frailty)
          frailty <- exp(2.5*fr$x)
    
    params <- list(
        n=2*n_frailty,
        S_ini=c(2*N*(1-f)*fr$p,2*N*f*fr$p),
        susceptibility=c(frailty, alpha*frailty))
     
  #print(params)
  if(!is.null(R)){
      beta <- get_beta(R, alpha, sd, f=f, N=N, n_frailty=n_frailty, gamma=gamma)
  }
    if(!is.null(k)){
      if(params$S_ini[abs(k)] <=1){
        return(rep(0, t))
      }
      if(k > 0){
          params$S_ini[k] <- params$S_ini[k] - 1
          params$S_ini[n_frailty + k] <- params$S_ini[n_frailty + k] + 1
      }else{
          params$S_ini[-k] <- params$S_ini[-k] + 1
          params$S_ini[n_frailty - k] <- params$S_ini[n_frailty - k] - 1
      
      }
      if(any(params$S_ini <0)){
          return(rep(0, t))
      }

    }
  #  print(params$S_ini)
        res <- run_det_cd(mm, rep(beta,t),params$S_ini + params$S_ini/sum(params$S_ini) , t, params$S_ini/sum(params$S_ini), susceptibility=params$susceptibility,gamma=gamma)
       
    if(!is.null(k)){
      if(k < 0){
          return(res$full_results[, paste("C[", - k, "]", sep="")]/params$S_ini[-k])
          
      }else{
          return(res$full_results[, paste("C[",n_frailty + k, "]", sep="")]/params$S_ini[n_frailty + k])
          
      
      }


    }
    }


    else{
        params <- list(
        n=2,
        S_ini=c(N, N),
        susceptibility=c(1, alpha))
        res <- run_det_cd(matrix(1, nrow=2, ncol=2)/2, rep(beta,t),params$S_ini , t, rep(1, 2), susceptibility=params$susceptibility,gamma=1/7)
    }

    
    
    exp <- res$full_results %>% as.data.frame() %>%select(starts_with("C"))
    if(sd>0){
        df <- data.frame(t=res$full_results[,1], vac=rowSums(exp[, (params$n/2+1):(params$n)]), unvac=rowSums(exp[, 1:(params$n/2)]))
        
    }else{
        df <- data.frame(t=res$full_results[,1], vac=exp[,2], unvac=exp[, 1])
    }
    N_vac <- 2*N*f
    N_unvac <- 2* N*(1-f)
    df$inc_unvac <- c(0, df$unvac[2:nrow(df)] - df$unvac[1:(nrow(df)-1)])
    df$inc_vac <- c(0, df$vac[2:nrow(df)] - df$vac[1:(nrow(df)-1)])
    df$HRR <- (df$inc_vac/(N_vac - df$vac))/(df$inc_unvac/(N_unvac - df$unvac))
    
    last_HRR_i <- max(which((N_vac-df$vac)> 1e-4 & (N_unvac-df$unvac)> 1e-4 & df$inc_unvac>1e-7 & df$inc_vac>1e-7))
    df$HRR[last_HRR_i:nrow(df)] <- df$HRR[last_HRR_i-1]
    return(list(sum=df %>% mutate(CRR=(vac/N_vac)/(unvac/N_unvac)), full=res$full_results))
}


get_beta <- function(R, alpha, sd, f=0.5, N=1000, n_frailty=100, gamma=1/2){
  fr <- get_frailty(sd=sd, n=n_frailty)
  frailty <- exp(2.5*fr$x)
  params <- list(
      n=2*n_frailty,
      S_ini=c(2*N*(1-f)*fr$p,2*N*f*fr$p),
      susceptibility=c(frailty, alpha*frailty))
  mm <- matrix(1, nrow=params$n, ncol=params$n)/params$n
  ng <- sweep(mm, MARGIN=1, params$susceptibility*params$S_ini/gamma/(2*N), `*`)
  eig <- Re(eigen(ng, only.values=T)$values[1])
  beta <- R / eig
  return(beta)
}

get_frailty <- function(mean=0.5, sd=1, n=100){
  variance <- sd^2
  common_factor <- (mean * (1 - mean) / variance) - 1
    shape1 <- mean * common_factor  # alpha
    shape2 <- (1 - mean) * common_factor  # beta
    

    #d <- distcrete::distcrete("beta",1, shape1=shape1, shape2=shape2)
    x <- (seq(0, 1, length.out=n+1) + 1/(n)/2)[1:n]

    p <- dbeta(x, shape1, shape2)
    return(list(x=x, p=p/sum(p)))
}





# Convert a binary contact matrix to a padded adjacency list suitable for
# det_model_adj.  Returns a list with:
#   neighbors  — integer matrix [max_degree x n], padded with 1
#   degree     — integer vector [n]
#   max_degree — scalar
contact_matrix_to_adj <- function(contact_matrix) {
  n <- nrow(contact_matrix)
  adj <- lapply(1:n, function(i) which(contact_matrix[i, ] != 0))
  degree <- lengths(adj)
  max_degree <- max(degree, 1L)
  neighbors <- matrix(1L,  nrow = max_degree, ncol = n)
  mask      <- matrix(0L,  nrow = max_degree, ncol = n)
  for (i in seq_len(n)) {
    if (degree[i] > 0) {
      neighbors[1:degree[i], i] <- adj[[i]]
      mask[1:degree[i], i]      <- 1L
    }
  }
  list(neighbors = neighbors, mask = mask, max_degree = max_degree)
}

run_det_cd <- function(mixing_matrix, beta_day, N, t, I_ini,
                      beta_norm=NULL,
                      susceptibility=NULL,
                      transmisibility=NULL, gamma=1/3,
                      waning=0,
                      import=0,
                      sparse=FALSE){
  if(is.null(beta_norm)) beta_norm <- N
  if(is.null(susceptibility)) susceptibility <- rep(1, dim(mixing_matrix)[1])
  if(is.null(transmisibility)) transmisibility <- rep(1, dim(mixing_matrix)[1])
  n <- dim(mixing_matrix)[1]

  if (sparse) {
    # odin2/dust2 path — O(edges) FOI via adjacency list
    adj <- contact_matrix_to_adj(mixing_matrix)
    # beta_day can be a vector; odin2 model takes a scalar beta so use mean
    beta_scalar <- mean(beta_day)
    params <- list(
      n=as.integer(n),
      max_degree=as.integer(adj$max_degree),
      beta=beta_scalar,
      gamma=gamma,
      waning=waning,
      neighbors=adj$neighbors,
      mask=adj$mask,
      susceptibility=susceptibility,
      transmisibility=transmisibility,
      S_ini=N - I_ini,
      I_ini=I_ini
    )
    sys <- dust2::dust_system_create(det_model_adj, params, time=0)
    dust2::dust_system_set_state_initial(sys)
    # state order: S[1..n], I[1..n], R[1..n], C[1..n]
    raw <- dust2::dust_system_simulate(sys, 1:t)  # [4n, t]
    S_mat <- raw[1:n, ]
    C_mat <- raw[(3*n+1):(4*n), ]
    # approximate n_SI from dC/dt (finite differences)
    nSI_mat <- cbind(C_mat[, 1, drop=FALSE],
                     C_mat[, 2:t] - C_mat[, 1:(t-1)])

    # build full_results data frame compatible with dense path
    full_res <- as.data.frame(t(rbind(S_mat, C_mat)))
    colnames(full_res) <- c(paste0("S[", 1:n, "]"), paste0("C[", 1:n, "]"))
    for (k in 1:n) full_res[[paste0("n_SI[", k, "]")]] <- nSI_mat[k, ]
    full_res[["t"]] <- 1:t

    C1   <- colSums(C_mat[1:(n/2), , drop=FALSE])
    C2   <- colSums(C_mat[(n/2+1):n, , drop=FALSE])
    S1   <- colSums(S_mat[1:(n/2), , drop=FALSE])
    S2   <- colSums(S_mat[(n/2+1):n, , drop=FALSE])
    nSI1 <- colSums(nSI_mat[1:(n/2), , drop=FALSE])
    nSI2 <- colSums(nSI_mat[(n/2+1):n, , drop=FALSE])

    main_comp <- data.frame(
      t        = 1:t,
      CRR      = (C2/sum(N[(n/2+1):n])) / (C1/sum(N[1:(n/2)])),
      HR       = (nSI2/S2) / (nSI1/S1),
      exposure = C1/sum(N[1:(n/2)]),
      unvac    = C1,
      vac      = C2
    )
    ind_comp <- list()
    for (k in 1:(n/2)) {
      tmp <- data.frame(
        t        = 1:t,
        exposure = C_mat[k, ] / N[k],
        CRR      = (C_mat[n/2+k, ]/N[n/2+k]) / (C_mat[k, ]/N[k]),
        HR       = (nSI_mat[n/2+k, ]/S_mat[n/2+k, ]) / (nSI_mat[k, ]/S_mat[k, ]),
        i        = as.character(k)
      )
      ind_comp[[k]] <- tmp
    }
    return(list(main=main_comp, ind=rbindlist(ind_comp, fill=TRUE),
                full_results=full_res))
  } else {
    params <- list(
      n=n,
      S_ini=N - I_ini,
      I_ini=I_ini,
      mixing_matrix=mixing_matrix,
      beta_day=beta_day,
      beta_norm=beta_norm,
      susceptibility=susceptibility,
      transmisibility=transmisibility,
      N_steps=t,
      waning=waning,
      interpolation_time=1:t,
      import=import,
      gamma=gamma
    )
    model <- det_model_cd$new(user=params)
  }
 
  res <- model$run(1:t)

  res <- as.data.frame(res)
  setDT(res)
  C1 <- rowSums(res[,paste("C[", 1:(n/2), "]", sep=""), with=FALSE])
  C2 <- rowSums(res[,paste("C[", (n/2+1):(n), "]", sep=""), with=FALSE])
  hazard1 <- rowSums(res[,paste("n_SI[", 1:(n/2), "]", sep=""), with=FALSE]) / rowSums(res[,paste("S[", 1:(n/2), "]", sep=""), with=FALSE])
  hazard2 <- rowSums(res[,paste("n_SI[", (n/2+1):(n), "]", sep=""), with=FALSE]) / rowSums(res[,paste("S[", (n/2+1):(n), "]", sep=""), with=FALSE])


  main_comp <- data.frame(t=1:t,
          CRR= (C2/sum(N[(n/2+1):n]))/(C1/sum(N[1:(n/2)])),
          HR=hazard2/hazard1,
          exposure=(C1/sum(N[1:(n/2)])),
          unvac=C1,
          vac=C2          
          )
          #,
          #type=rep(c("Control", "Vaccinated"), each=t))
ind_comp <- list()
for(i in 1:(n/2)){
    
    tmp <- data.frame(t=1:t,
      exposure=res[, paste("C[", i, "]", sep=""), with=FALSE]/N[i],
      CRR=(res[, paste("C[", (n/2+i), "]", sep=""), with=FALSE]/N[(n/2+i)])/(res[, paste("C[", i, "]", sep=""), with=FALSE]/N[i]),
      HR=(res[, paste("n_SI[", (n/2+i), "]", sep=""), with=FALSE]/res[, paste("S[", (n/2+i), "]", sep=""), with=FALSE])/(res[, paste("n_SI[", i, "]", sep=""), with=FALSE]/res[, paste("S[", i, "]", sep=""), with=FALSE]),
      i=as.character(i)) # %>%colnames() <- c("t", "CRR", "HR", "i")
      colnames(tmp) <- c("t","exposure", "CRR", "HR", "i")
  ind_comp[[i]] <- tmp  

  }
  return(list(main=main_comp,ind=rbindlist(ind_comp, fill=TRUE), full_results=as.data.frame(res)))
}



plot_cd <- function(res, var="I"){
  long <- res %>% select(t, starts_with(var)) %>% tidyr::pivot_longer(-t, names_to="group", values_to="value")

  ggplot(long) + geom_line(aes(x=t, y=value, color=group)) + ylab(var)
}




run_det_ncd <- function(N, t,  susceptibility){

  params <- list(
    n=length(N),
    S_ini=N,
    
    susceptibility=susceptibility
  )
  model <-det_model_ncd$new(user=params)
 
  res <- model$run(1:t)
 
  return(list(full_results=res))
}


#' Runs SIR model with given mixing matrix, beta over time and other key parameters
run_model <- function(mixing_matrix, beta_day, N, t, I_ini,
                      n_particles, n_threads=1, beta_norm=NULL,
                      susceptibility=NULL,
                      transmisibility=NULL, dt=0.2, print_index=FALSE, gamma=1/3,
                      import=0){
  if(is.null(beta_norm)) beta_norm <- N
  if(is.null(susceptibility)) susceptibility <- rep(1, dim(mixing_matrix)[1])
  if(is.null(transmisibility)) transmisibility <- rep(1, dim(mixing_matrix)[1])
  params <- list(
    n=dim(mixing_matrix)[1],
    S_ini=N - I_ini,
    I_ini=I_ini,
    mixing_matrix=mixing_matrix,
    beta_day=rep(beta_day, each=1/dt),
    beta_norm=beta_norm,
    susceptibility=susceptibility,
    transmisibility=transmisibility,
    N_steps=t/dt,
    dt=dt,
    import=import,
    gamma=gamma
  )
  dust_model <- odin_model$new(pars = params,
                               step = 1,
                               n_particles = n_particles,
                               n_threads = n_threads
                               )
  if(print_index){
    print(dust_model$info()$index)
  }
  res <- dust_model$simulate(1:t/dt)
  n <- dim(mixing_matrix)[1]
  R <- res[(2+2*n):(1+3*n),,t]
  fractions <- (R-I_ini) / (N + 1e-9)

  return(list(fractions=fractions,
              R=R,
              full_results=res))
}

#' Plot 2x2 matrix of relative risk
plot_2x2_RR <- function(fractions){

  RR <- fractions[2,]  / fractions[1,]
  ggplot(data.frame(RR=RR)) + geom_histogram(aes(x=RR))

}

#' Plot model trajectories
plot_history <- function(hist, var="I"){
  n <- (dim(hist)[1] - 1)/4
  print(n)
  sims <- dim(hist)[2]
  all_dfs <- list()
  for(i in 1:n){
    group <- hist[c((1 + i), (1+n+i), (1+2*n+i)),,]
    dfs <- list()
    for(sim in 1:sims){
      dfs[[sim]] <- data.frame(S=group[1, sim,],
                       I=group[2, sim,],
                       R=group[3, sim,],
                       t=1:dim(hist)[3],
                       sim=sim,
                       group=i)
    }
    
    all_dfs[[i]] <- rbindlist(dfs)
  }
  df <- rbindlist(all_dfs)
  df[, factor_group:=paste(sim, group)]
  
  ggplot(df) + geom_line(aes(x=t, y=get(var), group=factor_group, color=factor(group))) + ylab(var)

}

#' Construct the NGM matrix from a mixing matrix and other parameters
cij_NGM <- function(c_ij, N, susceptibility, transmisibility, gamma=1/3, norm_contacts=NULL){
  if(is.null(norm_contacts)){
    norm <- c_ij %*% N/sum(N)
  }else{
    N_conts <- as.numeric(norm_contacts %*% N)
    norm <- c_ij %*% N/sum(N)/N_conts*sum(N_conts)
  }

  c_ij <- c_ij/as.numeric(norm)
  NGM <- c_ij %*% diag(transmisibility)*N/sum(N)*susceptibility #NGM = suscept_i * S_i/N c_ij *trans_j
#  print(NGM)
  beta_R <- Re(eigen(NGM, only.values=T)$values[1]/gamma)
  return(list(c_ij=c_ij,
              NGM=NGM,
              beta_R=beta_R))
}

#' Run the model with 4 compartments for the paper
run_4x4 <- function(cimat, crmat, R0, susceptibility, transmisibility, gamma=1/3, n=100){
  N_N <- 90000
  N_I <- 10000
  N_IH <- 5000
  N_IL <- N_I - N_IH
  N_NH <- 0.1*N_N
  N_NL <- N_N - N_NH
  N <- c(N_NL, N_NH, N_IL, N_IH)  
  #cimat <- 
  cimat <- cimat/rowSums(cimat)
  #crmat <- matrix(c(1, 1, 1, 1), nrow=2)
  #crmat <- crmat/rowSums(crmat)
  
  c_ij <- matrix(1, nrow=4, ncol=4)
  c_ij[1:2, 1:2] <- cimat[1,1]*crmat
  c_ij[1:2, 3:4] <- cimat[1,2]*crmat
  c_ij[3:4, 1:2] <- cimat[2,1]*crmat
  c_ij[3:4, 3:4] <- cimat[2,2]*crmat

  norm_contacts <- cbind(rbind(crmat, crmat),rbind(crmat, crmat))

  input_mats <- cij_NGM(c_ij, N, susceptibility, transmisibility, norm_contacts=norm_contacts)
  beta <- R0/input_mats$beta_R
#  print("final c_ij")
#  print(input_mats$c_ij)
  res <- run_model(input_mats$c_ij, rep(beta, 200), N, 200, 0.001*N, n, 1,
                   susceptibility=susceptibility, transmisibility = transmisibility, gamma=gamma)
  res$N <- N
  return(res)

}

#' Run the model with 4 compartments and then either on the mean or for each simulation run a regression model to estimate regression coefficients
run_regs <- function(cimat, crmat, R0, susceptibility, transmisibility, on_mean=FALSE,
                     n=100, gamma=1/3){
  
  res <- run_4x4(cimat,
                 crmat,
                 R0,
                 susceptibility,
                 transmisibility,
                 n=n,
                 gamma=gamma)
  N <- res$N

  if(on_mean){
    data <- data.frame(N=N, I=rowMeans(res$full_results[10:13, , 200]), ethnicity=factor(c("N", "N", "I", "I"), levels=c("N", "I")), risk=c("L", "H", "L", "H"))
    return(glm(I~offset(log(N)) + ethnicity + risk, data=data, family=poisson))
  }else{
    sums <- list()
    for(i in 1:n){
      data <- data.frame(N=N, I=res$full_results[10:13,i, 200], ethnicity=c("N", "N", "I", "I"), risk=c("L", "H", "L", "H"))
      m <- glm(I~offset(log(N)) + ethnicity + risk, data=data, family=poisson)
      sums[[i]] <- m
    }
    return(sums)
  }
}


#' Theme for plots
add_theme <- function(q){
  q + theme_bw() + theme(text = element_text(size=8))+ scale_size_identity()
}



