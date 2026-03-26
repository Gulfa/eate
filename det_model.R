source("utils.R")

det_model_cd  <- odin::odin("det_mod_cd.R")
det_model_ncd <- odin::odin("det_mod_ncd.R")
det_model_adj <- odin2::odin("det_mod_adj.R")

# ---------------------------------------------------------------------------
# Core ODE runners
# ---------------------------------------------------------------------------

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
    vac=C2)

  ind_comp <- list()
  for(i in 1:(n/2)){
    tmp <- data.frame(t=1:t,
      exposure=res[, paste("C[", i, "]", sep=""), with=FALSE]/N[i],
      CRR=(res[, paste("C[", (n/2+i), "]", sep=""), with=FALSE]/N[(n/2+i)])/(res[, paste("C[", i, "]", sep=""), with=FALSE]/N[i]),
      HR=(res[, paste("n_SI[", (n/2+i), "]", sep=""), with=FALSE]/res[, paste("S[", (n/2+i), "]", sep=""), with=FALSE])/(res[, paste("n_SI[", i, "]", sep=""), with=FALSE]/res[, paste("S[", i, "]", sep=""), with=FALSE]),
      i=as.character(i))
    colnames(tmp) <- c("t","exposure", "CRR", "HR", "i")
    ind_comp[[i]] <- tmp
  }
  return(list(main=main_comp, ind=rbindlist(ind_comp, fill=TRUE), full_results=as.data.frame(res)))
}


run_det_ncd <- function(N, t, susceptibility){
  params <- list(
    n=length(N),
    S_ini=N,
    susceptibility=susceptibility
  )
  model <- det_model_ncd$new(user=params)
  res <- model$run(1:t)
  return(list(full_results=res))
}

# ---------------------------------------------------------------------------
# Frailty / heterogeneity model wrappers
# ---------------------------------------------------------------------------

run_frailty <- function(lambda, alpha, sd, N=1000, t=5000, n_frailty=200){
  if(sd>0){
    f <- get_frailty(sd=sd, n=n_frailty)
    frailty <- exp(log(lambda) + 2.5*f$x)
    params <- list(
      n=2*n_frailty,
      S_ini=N*rep(f$p, 2),
      susceptibility=c(frailty, alpha*frailty))
  }else{
    params <- list(
      n=2,
      S_ini=c(N, N),
      susceptibility=c(lambda, lambda*alpha))
  }
  model <- det_model_ncd$new(user=params)
  res <- model$run(0:t)
  exp <- res %>% as.data.frame() %>% select(starts_with("C"))
  if(sd>0){
    df <- data.frame(t=res[,1], vac=rowSums(exp[, (params$n/2+1):(params$n)]), unvac=rowSums(exp[, 1:(params$n/2)]))
  }else{
    df <- data.frame(t=res[,1], vac=exp[,2], unvac=exp[, 1])
  }
  df$inc_unvac <- rowSums(res[, paste("n_SR[", 1:(params$n/2), "]", sep="")])
  df$inc_vac   <- rowSums(res[, paste("n_SR[", (params$n/2+1):(params$n), "]", sep="")])
  df$HRR <- (df$inc_vac / (N - df$vac))/(df$inc_unvac / (N - df$unvac))

  strats <- list()
  for(i in 1:n_frailty){
    df_strat <- data.frame(t=res[,1], frailty_group=as.character(i))
    df_strat[,"vac"]      <- res[, paste("C[", n_frailty + i, "]", sep="")] / params$S_ini[n_frailty + i]
    df_strat[,"unvac"]    <- res[, paste("C[", i, "]", sep="")] / params$S_ini[i]
    df_strat$haz_unvac    <- res[, paste("hazzard_S[", i, "]", sep="")]
    df_strat$haz_vac      <- res[, paste("hazzard_S[", n_frailty + i, "]", sep="")]
    df_strat$HRR          <- (df_strat$haz_vac)/(df_strat$haz_unvac)
    df_strat$CRR          <- (df_strat$vac/df_strat$unvac)
    strats[[i]] <- df_strat
  }
  return(list(sum=df %>% mutate(CRR=vac/unvac, VE=1-CRR), strats=rbindlist(strats)))
}


run_frailty_cd <- function(alpha, sd, beta=1, R=NULL, f=0.5, N=1000, t=100, n_frailty=100, k=NULL, gamma=1/2){
  if(sd>0){
    fr <- get_frailty(sd=sd, n=n_frailty)
    frailty <- exp(2.5*fr$x)
    params <- list(
      n=2*n_frailty,
      S_ini=c(2*N*(1-f)*fr$p, 2*N*f*fr$p),
      susceptibility=c(frailty, alpha*frailty))

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
    res <- run_det_cd(mm, rep(beta,t), params$S_ini + params$S_ini/sum(params$S_ini), t, params$S_ini/sum(params$S_ini), susceptibility=params$susceptibility, gamma=gamma)

    if(!is.null(k)){
      if(k < 0){
        return(res$full_results[, paste("C[", -k, "]", sep="")]/params$S_ini[-k])
      }else{
        return(res$full_results[, paste("C[", n_frailty + k, "]", sep="")]/params$S_ini[n_frailty + k])
      }
    }
  }else{
    params <- list(
      n=2,
      S_ini=c(N, N),
      susceptibility=c(1, alpha))
    res <- run_det_cd(matrix(1, nrow=2, ncol=2)/2, rep(beta,t), params$S_ini, t, rep(1, 2), susceptibility=params$susceptibility, gamma=1/7)
  }

  exp <- res$full_results %>% as.data.frame() %>% select(starts_with("C"))
  if(sd>0){
    df <- data.frame(t=res$full_results[,1], vac=rowSums(exp[, (params$n/2+1):(params$n)]), unvac=rowSums(exp[, 1:(params$n/2)]))
  }else{
    df <- data.frame(t=res$full_results[,1], vac=exp[,2], unvac=exp[, 1])
  }
  N_vac <- 2*N*f
  N_unvac <- 2*N*(1-f)
  df$inc_unvac <- c(0, df$unvac[2:nrow(df)] - df$unvac[1:(nrow(df)-1)])
  df$inc_vac   <- c(0, df$vac[2:nrow(df)] - df$vac[1:(nrow(df)-1)])
  df$HRR <- (df$inc_vac/(N_vac - df$vac))/(df$inc_unvac/(N_unvac - df$unvac))
  last_HRR_i <- max(which((N_vac-df$vac)> 1e-4 & (N_unvac-df$unvac)> 1e-4 & df$inc_unvac>1e-7 & df$inc_vac>1e-7))
  df$HRR[last_HRR_i:nrow(df)] <- df$HRR[last_HRR_i-1]
  return(list(sum=df %>% mutate(CRR=(vac/N_vac)/(unvac/N_unvac)), full=res$full_results))
}

# ---------------------------------------------------------------------------
# Network / mean-field wrappers
# ---------------------------------------------------------------------------

run_mean_field <- function(beta=1, N=100, pl_alpha=3, alpha=1, t=100, vac_frac=0.5, vac=NULL, gamma=1/3, c_ij=NULL, k_mean=6){
  if(is.null(c_ij)){
    c_ij <- get_conact_matrix_pl(N, pl_alpha, mean_k=k_mean)
  }
  if(is.null(vac)){
    vac <- sample(1:N, vac_frac*N)
  }
  susept <- rep(1, N)
  susept[vac] <- alpha
  non_vac <- setdiff(1:N, vac)
  res <- run_det_cd(c_ij, rep(N*beta/k_mean, t), rep(1,N), t, I_ini=c(1,1,rep(0, N-2)), gamma=gamma, beta_norm=rep(1, N), susceptibility=susept, sparse=TRUE)

  sum <- data.frame(t=res$full_results[["t"]], vac=rowSums(res$full_results[, paste("C[", vac, "]", sep="")]), unvac=rowSums(res$full_results[, paste("C[", non_vac, "]", sep="")]))
  sum$CRR <- (sum$vac/(vac_frac*N))/(sum$unvac/((1-vac_frac)*N))
  return(list(sum=sum, full=res$full_results))
}

# ---------------------------------------------------------------------------
# EATE estimation
# ---------------------------------------------------------------------------

get_eate_frailty <- function(alpha, sd, beta=1, R=NULL, f=0.5, N=1000, t=100, n_frailty=100){
  full_res <- run_frailty_cd(alpha, sd, N=N, beta=beta, R=R, t=t, f=f, n_frailty=n_frailty)
  full_eff <- 0
  for(k in 1:n_frailty){
    res_k  <- run_frailty_cd(alpha, sd, N=N, beta=beta, R=R, t=t, f=f, n_frailty=n_frailty, k=k)
    res_mk <- run_frailty_cd(alpha, sd, N=N, beta=beta, R=R, t=t, f=f, n_frailty=n_frailty, k=-k)
    N_unvac <- full_res$full[, paste("S[", k, "]", sep="")][1]
    N_vac   <- full_res$full[, paste("S[", n_frailty + k, "]", sep="")][1]
    if(all(res_k==0) | all(res_mk==0)){
      next
    }
    eff  <- res_k / (full_res$full[ paste("C[", k, "]", sep="")]/N_unvac)
    eff2 <- (full_res$full[ paste("C[", n_frailty + k, "]", sep="")]/N_vac)/(res_mk)
    full_eff <- full_eff + eff*N_unvac + eff2 * N_vac
  }
  full_eff <- as.numeric(unlist(full_eff / (2*N)))
  return(full_res$sum %>% mutate(EATE=full_eff))
}


get_eate_network <- function(alpha=0.5, beta=1, R=NULL, f=0.5, N=200, t=15, pl_alpha=3, c_ij=NULL, n_vac=10){
  if(is.null(c_ij)){
    c_ij <- get_conact_matrix_pl(N, pl_alpha)
  }

  run_vac <- function(){
    vac <- sample(1:N, f*N)
    full_res <- run_mean_field(beta=beta, N=N, alpha=alpha, t=t, vac_frac=f, gamma=1, c_ij=c_ij, vac=vac)
    denom <- rep(0, t)
    num   <- rep(0, t)
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

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

plot_cd <- function(res, var="I"){
  long <- res %>% select(t, starts_with(var)) %>% tidyr::pivot_longer(-t, names_to="group", values_to="value")
  ggplot(long) + geom_line(aes(x=t, y=value, color=group)) + ylab(var)
}

# ---------------------------------------------------------------------------
# Analysis helpers (used in run_paper.R)
# ---------------------------------------------------------------------------

get_effects <- function(res){
  exp <- res %>% as.data.frame() %>% select(starts_with("C"))
  n_groups <- ncol(exp)
  Ns <- (res %>% as.data.frame() %>% select(starts_with("S")))[1,] %>% as.numeric()
  Ns <- Ns[seq(1, n_groups, by=2)]
  exp$unvac <- rowSums(exp[, seq(1, n_groups, by=2)])
  exp$vac   <- rowSums(exp[, seq(2, n_groups, by=2)])
  RR_i <- matrix(NA, nrow=nrow(exp), ncol=n_groups/2)
  for(i in 1:(n_groups/2)){
    RR_i[, i] <- exp[, paste0("C[", 2*i, "]")]/exp[, paste0("C[", 2*i-1, "]")]
  }
  AV_RR <- exp(Ns %*% t(log(RR_i))/sum(Ns)) %>% as.numeric()
  return(data.frame(t=res[, "t"], MarginalRatio=exp$vac/exp$unvac, AverageRatio=AV_RR))
}


plot_marg_vs_avg_4 <- function(lambda1, lambda2, alpha1, alpha2, N1, N2, t=150){
  params <- list(
    n=4,
    S_ini=c(N1/2, N1/2, N2/2, N2/2),
    susceptibility=c(lambda1, lambda1*alpha1, lambda2, lambda2*alpha2))
  model <- det_model_ncd$new(user=params)
  res <- model$run(0:t)
  df <- get_effects(res)
  exp <- res %>% as.data.frame() %>% select("t", starts_with("C")) %>%
    tidyr::pivot_longer(cols=-t) %>%
    mutate(vaccinated=grepl("C\\[2\\]|C\\[4\\]", name),
           name=recode(name, `C[1]`="Unvac, high risk", `C[2]`="Vac, high risk",
                       `C[3]`="Unvac, low risk", `C[4]`="Vac, low risk"))
  exp$value <- exp$value / rep(params$S_ini, nrow(exp)/length(params$S_ini))
  q1 <- ggplot(exp) + geom_line(aes(x=t, y=value, col=name), size=1.7) + theme_minimal() + labs(y="Attack rate", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("", palette="Dark2") + theme(legend.position="bottom")
  q2 <- ggplot(df %>% tidyr::pivot_longer(cols=-t)) + geom_line(aes(x=t, y=1-value, col=name), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Effect measure", palette="Dark2") + theme(legend.position="bottom") + scale_y_continuous(labels=scales::percent_format(accuracy=1))
  main <- cowplot::plot_grid(q1, q2, ncol=2, labels=c("A", "B"), rel_widths=c(1,1))
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      glue::glue("Lambda1={lambda1}, Lambda2={lambda2}, Alpha1={alpha1}, Alpha2={alpha2}, N1={N1}, N2={N2}"),
      fontface='bold', x=0, hjust=0) +
    theme(plot.margin=margin(0, 0, 0, 7))
  cowplot::plot_grid(title, main, ncol=1, rel_heights=c(0.1, 1))
}


run_nl <- function(frac_vac=0.5, N=1000, beta=3, alpha=0.3){
  cd <- run_det_cd(mixing_matrix=matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2), beta_day=rep(beta, 50),
                   t=50, N=c(N*(1-frac_vac), N*frac_vac), I_ini=c(1,1),
                   susceptibility=c(1, alpha), gamma=1)$full_results
  out <- data.frame(t=rep(cd[, "t"], 2),
                    a=c(cd[, "C[1]"]/(N*(1-frac_vac)), cd[, "C[2]"]/(N*frac_vac)),
                    group=rep(c("Control", "Vaccinated"), each=length(cd[,"t"])))
  return(out)
}


get_average_observed <- function(frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1)){
  ao <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-frac), frac*N), 500, c(0.1,0), susceptibility=susceptibility, gamma=gamma)
  return(data.frame(t=ao$full_results[, "t"],
                    IR=ao$full_results[, "C[2]"]/((frac)*N) / (ao$full_results[, "C[1]"]/((1-frac)*N)),
                    attack_rate=ao$full_results[, "C[1]"]/((1-frac)*N)))
}


get_individual_effect <- function(vac_frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1)){
  full_1 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac)-1, N*vac_frac+1), 500, c(0.1,0), susceptibility=susceptibility, gamma=gamma)
  full_0 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac),   N*vac_frac),   500, c(0.1,0), susceptibility=susceptibility, gamma=gamma)
  return(data.frame(t=full_1$main[, "t"],
                    attack_rate=full_0$main$unvac/(N*(1-vac_frac)),
                    IR=full_1$full_results[, "C[2]"]/(N*vac_frac+1)/(full_0$full_results[, "C[1]"]/(N*(1-vac_frac)))))
}


get_EATE <- function(vac_frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1), t=500){
  full_1  <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,t), c(N*(1-vac_frac)-1, N*vac_frac+1), t, c(0.1,0), susceptibility=susceptibility, gamma=gamma)
  full_0  <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,t), c(N*(1-vac_frac),   N*vac_frac),   t, c(0.1,0), susceptibility=susceptibility, gamma=gamma)
  full_m1 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,t), c(N*(1-vac_frac)+1, N*vac_frac-1), t, c(0.1,0), susceptibility=susceptibility, gamma=gamma)

  num   <- (1-vac_frac)*full_1$full_results[, "C[2]"] + vac_frac*full_0$full_results[, "C[2]"]
  denom <- (1-vac_frac)*full_0$full_results[, "C[1]"] + vac_frac*full_m1$full_results[, "C[1]"]

  return(data.frame(t=full_1$main[, "t"],
                    attack_rate=full_0$main$unvac/(N*(1-vac_frac)),
                    eate=num/denom,
                    CRR=full_0$full_results[, "C[2]"]/(N*vac_frac) / (full_0$full_results[, "C[1]"]/(N*(1-vac_frac)))))
}


get_EATE_frailty <- function()


get_HH <- function(vac_frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1)){
  full_1 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac)-1, N*vac_frac+1), 500, c(0.1,0), susceptibility=susceptibility, gamma=gamma)
  full_0 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac),   N*vac_frac),   500, c(0.1,0), susceptibility=susceptibility, gamma=gamma)
  return(data.frame(t=full_1$main[, "t"],
                    attack_rate=full_0$main$unvac/(N*(1-vac_frac)),
                    IR=full_1$full_results[, "C[2]"]/(N*vac_frac+1)/(full_0$full_results[, "C[1]"]/(N*(1-vac_frac)))))
}


def_run_4_4_nl <- function(f_c, f_v, lambda_h, lambda_l, alpha, beta=0.3, N=10000, model=NULL, M=NULL){
  t <- 200
  res <- run_det_cd(M, rep(beta, t), N=c(N*f_c, N*(1-f_c), N*f_v, N*(1-f_v)), t, c(1,1,1,1),
                   susceptibility=c(lambda_h, lambda_l, lambda_h*alpha, lambda_l*alpha), gamma=1/7)
  res <- res$full_results
  df <- data.frame(t=res[,1],
                   CRR_t=(N-2 - res[,4]-res[,5])/(N-2 - res[,2]-res[,3]),
                   CRR_h=((N*f_v-1 - res[,4])/(N*f_v-1))/((N*f_c-1-res[,2])/(N*f_c-1)),
                   CRR_l=((N*(1-f_v)-1 - res[,5])/(N*(1-f_v)-1))/((N*(1-f_c)-1-res[,3])/(N*(1-f_c)-1))) %>% filter(t>1)
  return(list(res=res, rr=df, N=c(N*f_c, N*(1-f_c), N*f_v, N*(1-f_v))))
}


strat_est <- function(res, frac_h=0.5){
  return(data.frame(t=res$rr$t, rr=res$rr$CRR_h*frac_h + res$rr$CRR_l*(1-frac_h)))
}
