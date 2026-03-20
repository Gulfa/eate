library(dplyr)
library(data.table)
library(ggplot2)
library(adaptivetau)

source("model.R")

#FIG 1


lambda1 <- 0.1
lambda2 <- 0.05
alpha <- 0.5
t <- seq(0, 200, by=0.05)

IRR1 = (1 - exp(-alpha*lambda1*t)) /(1 - exp(-lambda1*t) )
IRR2 = (1 - exp(-alpha*lambda2*t)) /(1 - exp(-lambda2*t) )

ggplot(data.frame(t=rep(t,2), IRR=c(IRR1, IRR2), lambda=rep(c(lambda1, lambda2), each=length(t))), aes(x=t, y=IRR, color=factor(lambda))) + geom_line(size=2) + theme_minimal() + theme(legend.position="top") + labs(x="Time", y="Incidence Rate Ratio") + scale_color_brewer("Lambda", palette="Dark2") + theme(text = element_text(size=24))

ggsave("article/figures/IRR.png", width=10, height=6, bg="white")


lambda <- 0.05
alpha <- 0.5
N <- 1000

frailties <- list()
for(lambda in c(0.05, 0.1)){
    for(sd in c(0.01, 0.1, 1/sqrt(12), 0.35, 0.4)){
        frailties[[length(frailties) + 1]] <- run_frailty(lambda, alpha, sd) %>% mutate(sd=round(sd,2), alpha=alpha, lambda=lambda)
    }
}
frailties <- rbindlist(frailties)

frailties$lambda <- factor(frailties$lambda, levels=c(0.05, 0.1), labels=c("Lambda=0.05", "Lambda=0.1"))

q3 <- ggplot(frailties) + geom_line(aes(x=unvac/1e3, y=vac/1e3, group = paste(lambda, sd), col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Cumulative incidence vaccinated", x="Cumulative incidence control") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")
q2 <- ggplot(frailties %>% filter(t< 150)) + geom_line(aes(x=t, y=HRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Hazard rate ratio", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + facet_grid(.~lambda)
q1 <- ggplot(frailties %>% filter(t< 150)) + geom_line(aes(x=t, y=CRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Cumulative incidence rate ratio", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + facet_grid(.~lambda)

cowplot::plot_grid(q1, q2,q3, ncol=1, labels=c("A", "B", "C"), rel_heights = c(1,1,1))

ggsave("output/attack_rate_hetero.png", width=12, height=16)


get_effects <- function(res){
    exp <- res %>% as.data.frame() %>%select(starts_with("C"))
    
    n_groups <- ncol(exp)
    Ns <- (res %>% as.data.frame() %>%select(starts_with("S")))[1,] %>% as.numeric()
    Ns <- Ns[seq(1, n_groups, by=2)] 
    exp$unvac <- rowSums(exp[, seq(1, n_groups, by=2)])
    exp$vac <- rowSums(exp[, seq(2, n_groups, by=2)])
    RR_i <- matrix(NA, nrow=nrow(exp), ncol=n_groups/2)
    for(i in 1:(n_groups/2)){
        RR_i[, i] <- exp[, paste0("C[", 2*i, "]")]/exp[, paste0("C[", 2*i-1, "]")]
    }
    AV_RR <- exp(Ns %*% t(log(RR_i))/sum(Ns)) %>% as.numeric()
    return(data.frame(t=res[, "t"], MarginalRatio=exp$vac/exp$unvac,
                                AverageRatio = AV_RR))

}


plot_marg_vs_avg_4 <- function(lambda1, lambda2, alpha1, alpha2, N1, N2, t=150) {
    params <- list(

        n=4,
        S_ini=c(N1/2, N1/2, N2/2, N2/2),
        susceptibility=c(lambda1, lambda1*alpha1, lambda2, lambda2*alpha2))
    model <-det_model_ncd$new(user=params)
    res <- model$run(0:t)
    df <- get_effects(res)
    exp <- res %>% as.data.frame() %>%select("t", starts_with("C"))  %>% tidyr::pivot_longer(cols=-t) %>% mutate(vaccinated=grepl("C\\[2\\]|C\\[4\\]", name), name=recode(name, `C[1]`="Unvac, high risk", `C[2]`="Vac, high risk", `C[3]`="Unvac, low risk", `C[4]`="Vac, low risk"))
    #calculate factional attack rate for each group dividing by params$S_ini
    exp$value <- exp$value / rep(params$S_ini, nrow(exp)/length(params$S_ini))
    q1 <- ggplot(exp) + geom_line(aes(x=t, y=value, col=name), size=1.7) + theme_minimal() + labs(y="Attack rate", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("", palette="Dark2") + theme(legend.position="bottom")
    q2 <- ggplot(df %>% tidyr::pivot_longer(cols=-t)) + geom_line(aes(x=t, y=1 - value, col=name), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Effect measure", palette="Dark2") + theme(legend.position="bottom")  + scale_y_continuous(labels=scales::percent_format(accuracy=1))
    main <- cowplot::plot_grid(q1, q2, ncol=2, labels=c("A", "B"), rel_widths = c(1,1))

    title <- cowplot::ggdraw() + 
  cowplot::draw_label(
   glue::glue("Lambda1={lambda1}, Lambda2={lambda2}, Alpha1={alpha1}, Alpha2={alpha2}, N1={N1}, N2={N2}"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
cowplot::plot_grid(
  title, main,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)
}

plot_marg_vs_avg_4(0.05, 0.02, 0.5, 0.5, 1000, 1000, t=300)
ggsave("output/marginal_vs_average_4groups_1.png", width=12, height=6)
plot_marg_vs_avg_4(0.05, 0.01, 10, 0.5, 20, 1000, t=300)
ggsave("output/marginal_vs_average_4groups_2.png", width=12, height=6)



res <- run_frailty(0.01, 0.5, 0.4, N=1000, t=250, n_frailty=20)

ggplot(res$strats) + geom_line(aes(x=unvac, y=vac, color=frailty_group))

ggplot(res$strats) + geom_line(aes(x=unvac, y=vac/unvac, color=frailty_group)) + geom_line(aes(x=unvac/1000, y=vac/unvac), data=res$sum) + theme_minimal() + labs(y="Cumulative incidence ratio", x="Cumulative incidence control") + theme(text = element_text(size=20)) + scale_color_brewer("Frailty group", palette="Dark2") + theme(legend.position="bottom")

ggplot(res$strats) + geom_line(aes(x=t, y=HRR, color=frailty_group)) + geom_line(aes(x=t, y=HRR), data=res$sum)


res2 <- run_frailty(0.01, 0.5, 0.15, N=1000, t=250, n_frailty=20)
res3 <- run_frailty(0.01, 0.5, 0.25, N=1000, t=250, n_frailty=20)
res4 <- run_frailty(0.01, 0.5, 0.0001, N=1000, t=250, n_frailty=3)
comb <- rbind(res$sum %>% mutate(sd=0.4), res2$sum %>% mutate(sd=0.15), res3$sum %>% mutate(sd=0.25), res4$sum %>% mutate(sd=0.00001))
ggplot(comb) + geom_line(aes(x=t, y=HRR, color=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Marginal Hazard rate ratio", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + theme(legend.position="bottom")
ggsave("output/hazard_rate_hetero.png", width=10, height=6)
ggplot(comb) + geom_line(aes(x=unvac/1000, y=CRR, color=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Marginal Cumulative incidence rate ratio", x="Fraction infected in control") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + theme(legend.position="bottom")
ggsave("output/cumulative_incidence_hetero.png", width=10, height=6)
# 6 groups linear model: 
lambda1 <- 0.003
lambda2 <- 0.001
lambda3 <- 0.05
alpha <- 0.5
params <- list(

        n=6,
        S_ini=c(N/2, N/2, N/2, N/2, N/2, N/2),
        susceptibility=c(lambda1, lambda1*alpha, lambda2, lambda2*alpha, lambda3, lambda3*alpha))
model <-det_model_ncd$new(user=params)
res <- model$run(0:500)
exp <- res %>% as.data.frame() %>%select(starts_with("C"))
exp$unvac <- exp$`C[1]` + exp$`C[3]` + exp$`C[5]`
exp$vac <- exp$`C[2]` + exp$`C[4]` + exp$`C[6]`

exp$inc_unvac <- c(0, exp$unvac[2:nrow(exp)] - exp$unvac[1:(nrow(exp)-1)])
exp$inc_vac <- c(0, exp$vac[2:nrow(exp)] - exp$vac[1:(nrow(exp)-1)])

   

df <- data.frame(t=res[, "t"], MarginalRatio=(exp[, "C[2]"] + exp[, "C[4]"] / (exp[, "C[1]"] + exp[, "C[3]"])),# + exp[, "C[3]"] + exp[, "C[5]"]),
                                AverageRatio = (exp[, "C[2]"] / exp[, "C[1]"] + exp[, "C[4]"] / exp[, "C[3]"]))# + exp[, "C[6]"] / exp[, "C[5]"])/3,
 #                               CumDiff = (exp[, "C[1]"] + exp[, "C[3]"] + exp[, "C[5]"] - exp[, "C[2]"] - exp[, "C[4]"] - exp[, "C[6]"]),
  #                              AverageDiff = (exp[, "C[1]"] - exp[, "C[2]"] + exp[, "C[3]"] - exp[, "C[4]"] + exp[, "C[5]"] - exp[, "C[6]"]),
                             #   MarginalHazardRatio=(exp$inc_vac/(3*N/2 - exp$vac))/(exp$inc_unvac/(3*N/2 - exp$unvac))) 


ggplot(df %>% tidyr::pivot_longer(cols=-t)) %>% mutate(nice=recode(name, CumIncRatio="Cumulative Incidence Ratio", HazardRateRatio="Hazard Rate Ratio"))) + geom_line(aes(x=t, y=1- value, col=nice), size=1.7) + theme_minimal() + labs(y="Vaccine Effectiveness(VE)", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Effect measure", palette="Dark2") + theme(legend.position="bottom") + scale_y_continuous(labels=scales::percent_format(accuracy=1))

ggsave("output/average_vs_overall.png", width=10, height=6)



res <- list()

for(sd in c(0.01, 0.4)){
    beta <- get_beta(2.5, 1, sd=sd, f=0.5, N=1000, n_frailty=25, gamma=1/2)
    print(beta)
    for(f in c(0.1, 0.7)){
        res[[length(res) + 1]] <- get_eate_frailty( 0.5,sd=sd,f=f,beta=beta, N=1000,t=100, n_frailty=25) %>% mutate(sd=sd, f=glue::glue("{round(f*100,0)}% vaccinated"))
    }
    
}
res <- rbindlist(res)
    


get_eate_frailty( 0.5,sd=sd,f=0.7,beta=beta, N=1000,t=100, n_frailty=25) %>% mutate(sd=sd, f=glue::glue("{round(f*100,0)}% vaccinated"))

q1 <- ggplot(res %>% filter(t<=100)) + geom_line(aes(x=t, y=EATE, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="EATE", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + facet_wrap(.~f)
q2 <- ggplot(res %>% filter(t<=100)) + geom_line(aes(x=t, y=CRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="CIR", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")+ facet_wrap(.~f)
q3 <- ggplot(res %>% filter(t<=100)) + geom_line(aes(x=t, y=HRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="HRR", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")+ facet_wrap(.~f)
cowplot::plot_grid(q1, q2, q3, ncol=1, labels=c("A", "B", "C"), rel_heights = c(1,1,1))

ggsave("output/effect_sir.png", width=10, height=12)


run_nl <- function(frac_vac=0.5, N=1000,beta=3, alpha=0.3){


    cd <- run_det_cd(mixing_matrix=matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2), beta_day=rep(beta, 50), t=50, N=c(N*(1-frac_vac),N*frac_vac), I_ini=c(1,1), susceptibility=c(1, alpha), gamma=1)$full_results
    out <- data.frame(t=rep(cd[, "t"], 2),
            a=c(cd[, "C[1]"]/(N*(1-frac_vac)), cd[, "C[2]"]/(N*frac_vac)), group=rep(c("Control", "Vaccinated"), each=length(cd[,"t"])))
    return(out)
}

run_nl(beta=10)

df <- rbind(run_nl(0.25) %>% mutate(frac_vac=glue::glue("25% vaccinated")), run_nl(0.01) %>% mutate(frac_vac=glue::glue("1% vaccinated")), run_nl(0.5) %>% mutate(frac_vac=glue::glue("50% vaccinated")))


df_eate <- rbind(
    get_EATE(0.01, 3) %>% mutate(frac_vac=glue::glue("1% vaccinated")),
    get_EATE(0.25, 3) %>% mutate(frac_vac=glue::glue("25% vaccinated")),
    get_EATE(0.5, 3) %>% mutate(frac_vac=glue::glue("50% vaccinated"))
)

p1 <- ggplot(df) + geom_line(aes(x=t, y=a, col=group), size=1.7) + theme_minimal() + labs(y="Attack rate", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Group", palette="Dark2") + facet_wrap(~frac_vac)

p2 <- ggplot(df_eate %>% filter(t<=50)) + geom_line(aes(x=t, y=IR, col=frac_vac), size=1.7) + theme_minimal() + labs(y="EATE", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Vaccination Coverage", palette="Set1")

cowplot::plot_grid(p1, p2, ncol=1, labels=c("A", "B"), rel_heights = c(1,0.5))

ggsave("output/attack_rate_nl.png", width=10, height=6)



c_ij <- get_conact_matrix_pl(300, 1.5, mean_k=6)


random_mixing  <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(2, 20), c(100, 100), 20, c(2,0), susceptibility=c(1,0.2), gamma=1)

CMMs <- list()
for(i in 1:20){
    print(i)
    r <- run_mean_field(beta=2, N=300, alpha=0.2, t=20, gamma=1, c_ij=c_ij)
    CMMs[[length(CMMs) + 1]] <- r$sum %>% select(t, CRR, vac, unvac) %>% mutate(run=i)
}
res <- rbindlist(CMMs)
ggplot(res) + geom_line(aes(x=t, y=CRR, group=run), alpha=0.5) + theme_minimal() + labs(y="Cumulative incidence rate ratio", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Run", palette="Dark2") + theme(legend.position="bottom") + geom_line(aes(x=t, y=CRR), color="red", size=1.7, data=random_mixing$main) + geom_line(aes(x=t, y=CRR), color="blue", size=1.7, data=res %>% group_by(t) %>% summarise(CRR=mean(CRR)))


ggplot(res) + geom_line(aes(x=t, y=vac + unvac, group=run), alpha=0.5) + theme_minimal() + labs(y="Cumulative incidence", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Run", palette="Dark2") + theme(legend.position="bottom") + geom_line(aes(x=t, y=vac+unvac), color="red", size=1.7, data=random_mixing$main)



get_average_observed <- function(frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1)){
    ao <-  run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-frac), frac*N), 500, c(0.1,0), susceptibility=susceptibility,gamma=gamma)
    return(data.frame(t=ao$full_results[, "t"], IR=ao$full_results[, "C[2]"]/((frac)*N) /( ao$full_results[, "C[1]"]/((1-frac)*N)), attack_rate=ao$full_results[, "C[1]"]/((1-frac)*N)))
}

get_individual_effect <- function(vac_frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1)){

    full_1 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac) - 1,N*vac_frac + 1), 500, c(0.1,0), susceptibility=susceptibility,gamma=gamma)
    full_0 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac) ,N*vac_frac), 500, c(0.1,0), susceptibility=susceptibility,gamma=gamma)

    return(data.frame(t=full_1$main[, "t"], attack_rate = full_0$main$unvac/(N*(1-vac_frac)), IR=full_1$full_results[, "C[2]"]/(N*vac_frac + 1)/(full_0$full_results[, "C[1]"]/(N*(1-vac_frac)))))
}


get_EATE <- function(vac_frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1), t=500){

    full_1 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,t), c(N*(1-vac_frac) - 1,N*vac_frac + 1), t, c(0.1,0), susceptibility=susceptibility,gamma=gamma)
    full_0 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,t), c(N*(1-vac_frac) ,N*vac_frac), t, c(0.1,0), susceptibility=susceptibility,gamma=gamma)
    full_m1 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,t), c(N*(1-vac_frac)+1 ,N*vac_frac-1), t, c(0.1,0), susceptibility=susceptibility,gamma=gamma)

    num <- (1-vac_frac)*full_1$full_results[, "C[2]"] + vac_frac*full_0$full_results[, "C[2]"]
    denom <- (1-vac_frac)*full_0$full_results[, "C[1]"] + vac_frac*full_m1$full_results[, "C[1]"]
    


    return(data.frame(t=full_1$main[, "t"], attack_rate = full_0$main$unvac/(N*(1-vac_frac)), eate=num/denom, CRR=full_0$full_results[, "C[2]"]/(N*vac_frac )/(full_0$full_results[, "C[1]"]/(N*(1-vac_frac)))))
}

get_EATE_frailty <- function()



get_HH <- function(vac_frac, beta, gamma=0.2, N=1000, susceptibility=c(1,0.1)){
    full_1 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac) - 1,N*vac_frac + 1), 500, c(0.1,0), susceptibility=susceptibility,gamma=gamma)
    full_0 <- run_det_cd(matrix(1, nrow=2, ncol=2), rep(beta,500), c(N*(1-vac_frac) ,N*vac_frac), 500, c(0.1,0), susceptibility=susceptibility,gamma=gamma)
    return(data.frame(t=full_1$main[, "t"], attack_rate = full_0$main$unvac/(N*(1-vac_frac)), IR=full_1$full_results[, "C[2]"]/(N*vac_frac + 1)/(full_0$full_results[, "C[1]"]/(N*(1-vac_frac)))))
}

get_indivudal_effect(0, 0.4)


a <-   get_individual_effect(0.5, 0.4, N=10)
b <-   get_EATE(0.5, 2, N=100, gamma=1, susceptibility=c(1,0.5), t=25)


c <- get_eate_network(alpha=0.5, beta=2, N=100, f=0.5,t=15, pl_alpha = 1, n_vac=30)

ggplot(c %>% tidyr::pivot_longer(cols=c(eate, CRR), names_to="measure", values_to="value")) + geom_line(aes(x=t, y=value, group=sim), size=1.7) + theme_minimal() + labs(y="Effect measure", x="Time") + theme(text = element_text(size=20))  + theme(legend.position="bottom") + facet_wrap(~measure)

ggplot(c) + geom_line(aes(x=t, y=eate, group=sim), size=1.7) + theme_minimal() + labs(y="EATE", x="Time") + theme(text = element_text(size=20)) + geom_line(aes(x=t, y=CRR, group=sim), color="red", size=1.7)  + geom_line(aes(x=t, y=eate), color="blue", size=1.7, data=b) 


ggplot(c) + geom_line(aes(x=t, y=eate/CRR, group=sim), size=1.7) + theme_minimal() + labs(y="EATE/CRR", x="Time") + theme(text = element_text(size=20))   + geom_line(aes(x=t, y=eate/CRR), color="blue", size=1.7, data=b) 


exp <- c %>% group_by(t) %>% summarise(eate=mean(eate), CRR=mean(CRR))

ggplot(exp) + geom_line(aes(x=t, y=eate/CRR), size=1.7) + theme_minimal() + labs(y="EATE/CRR", x="Time") + theme(text = element_text(size=20))   + geom_line(aes(x=t, y=eate/CRR), color="blue", size=1.7, data=b)

a$IR/b$IR

res <-rbind(
            get_average_observed(0.2, 0.4) %>% mutate(eff="Observed, coverage=20%", beta="beta=2", N=1000, cov=0.2, type="Observed"),
            get_average_observed(0.5, 0.4) %>% mutate(eff= "Observed, coverage=50%", beta="beta=2", N=1000,cov=0.5, type="Observed"),
            get_EATE(0.2, 0.4) %>% mutate(eff="EATE coverage=20%", beta="beta=2", N=1000,cov=0.2, type="EATE"),
            get_EATE(0.5, 0.4) %>% mutate(eff= "EATE coverage=50%", beta="beta=2", N=1000, cov=0.5, type="EATE"),
            get_average_observed(0.2, 0.6) %>% mutate(eff="Observed, coverage=20%", beta="beta=3", N=1000, cov=0.2, type="Observed"),
            get_average_observed(0.5, 0.6) %>% mutate(eff= "Observed, coverage=50%", beta="beta=3", N=1000, cov=0.5, type="Observed"),
            get_EATE(0.2, 0.6) %>% mutate(eff="EATE coverage=20%", beta="beta=3", N=1000, cov=0.2, type="EATE"),
            get_EATE(0.5, 0.6) %>% mutate(eff= "EATE coverage=50%", beta="beta=3", N=1000, cov=0.5, type="EATE"),
            get_average_observed(0.2, 0.4, N=10) %>% mutate(eff="Observed, coverage=20%", beta="beta=2", N=10, cov=0.2, type="Observed"),
            get_average_observed(0.5, 0.4, N=10) %>% mutate(eff= "Observed, coverage=50%", beta="beta=2", N=10, cov=0.5, type="Observed"),
            get_EATE(0.2, 0.4, N=10) %>% mutate(eff="EATE coverage=20%", beta="beta=2", N=10, cov=0.2, type="EATE"),
            get_EATE(0.5, 0.4, N=10) %>% mutate(eff= "EATE coverage=50%", beta="beta=2", N=10,     cov=0.5, type="EATE"),
            get_average_observed(0.2, 0.6, N=10) %>% mutate(eff="Observed, coverage=20%", beta="beta=3", N=10, cov=0.2, type="Observed"),
            get_average_observed(0.5, 0.6, N=10) %>% mutate(eff= "Observed, coverage=50%", beta="beta=3", N=10, cov=0.5, type="Observed"),
            get_EATE(0.2, 0.6, N=10) %>% mutate(eff="EATE coverage=20%", beta="beta=3", N=10, cov=0.2, type="EATE"),
            get_EATE(0.5, 0.6, N=10) %>% mutate(eff= "EATE coverage=50%", beta="beta=3", N=10, cov=0.5, type="EATE")
)
            
res$N <- paste("N =", res$N)

ggplot(res %>% filter(t<100)) + geom_line(aes(x=t, y=IR, col=eff), size=1.7) + theme_minimal() + labs(y="Effect measure", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Effect", palette="Dark2") + facet_grid(N~beta)
ggsave("output/different_effects.png", width=10, height=6)

ggplot(res %>% filter(t<200 & grepl("EATE", eff))) + geom_line(aes(x=attack_rate, y=IR, col=eff), size=1.7) + theme_minimal() + labs(y="Effect measure", x="Cummulative Incidence in control group") + theme(text = element_text(size=20)) + scale_color_brewer("Effect", palette="Dark2") + facet_grid(N~beta)
ggsave("output/different_effect_atk.png", width=10, height=6)


ggplot(res %>%filter(type=="EATE" & N==1000 & t<=200), aes(x=t, y=IR, col=factor(cov))) + geom_line(size=1.7) + theme_minimal() + labs(y="EATE(rate)", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Coverage", palette="Dark2") + facet_grid(.~beta)
ggsave("output/different_effects_cov.png", width=10, height=6)



a1 <-  get_EATE(0.7, 0.6, N=10)
a2 <-  get_EATE(0.7, 0.6, N=1000) 
b1 <- get_average_observed(0.7, 0.6, N=10) 
b2 <- get_average_observed(0.7, 0.6)


a1 <-  get_EATE(0.1, 0.6, N=10)
a2 <-  get_EATE(0.1, 0.6, N=1000) 
b1 <- get_average_observed(0.1, 0.6, N=10) 
b2 <- get_average_observed(0.1, 0.6)

overall <- list()
for( vac_frac in c(0.2, 0.5, 0.7)){
    eate <- c()
    obs <- c()
    for(N in c(10, 1000)){
        a1 <-  get_EATE(vac_frac, 0.6, N=N)
        b1 <- get_average_observed(vac_frac, 0.6, N=N) 
        eate <- c(eate, a1$IR)
        obs <- c(obs, b1$IR)
    }
    overall[[length(overall) + 1]] <- data.frame(eate=eate, obs=obs, N=rep(c(10, 1000), each=length(eate)/2), vac_frac=vac_frac)
}


overall <- rbindlist(overall)


ggplot(overall,aes(x=obs, y=eate, col=factor(N))) + geom_point(size=3) + theme_minimal() + labs(y="EATE", x="Observed difference (ATE)") + theme(text = element_text(size=20)) + scale_color_brewer("N", palette="Dark2") + geom_abline(slope=1, intercept=0, linetype="dashed") + facet_wrap(.~vac_frac, scales="free_x")

ggsave("output/eate_obs.png", width=12, height=8)



def_run_4_4_nl<- function(f_c, f_v, lambda_h, lambda_l, alpha, beta=0.3, N=10000, model=det_model_non_linear, M=NULL){

    t <- 200
    res <- run_det_cd(M,rep(beta, t), N=c(N*f_c, N*(1-f_c), N*f_v, N*(1-f_v)), t, c(1,1,1,1), susceptibility=c(lambda_h, lambda_l, lambda_h*alpha, lambda_l*alpha),gamma=1/7)

    
    res <- res$full_results
    df <- data.frame(t=res[,1], CRR_t=(N-2 - res[,4]-res[,5])/(N-2 - res[,2]-res[,3]), CRR_h=((N*f_v-1 - res[,4])/(N*f_v-1))/((N *f_c-1- res[,2])/(N*f_c-1)),CRR_l=((N*(1-f_v)-1 - res[,5])/(N*(1-f_v)-1))/((N *(1-f_c)-1- res[,3])/(N*(1-f_c)-1))) %>% filter(t>1)
    
    return(list(res=res, rr=df, N=c(N*f_c, N*(1-f_c), N*f_v, N*(1-f_v))) )
}

lambda_h <- 3
lambda_l <- 1
alpha <- 0.5


m_random <- matrix(rep(1, 16), nrow=4, ncol=4)
m_sep <- matrix(c(1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1), nrow=4, ncol=4, byrow=TRUE)

non <- def_run_4_4_nl(0.5, 0.5, lambda_h, lambda_l, alpha, model=det_non_linear, M=m_random)
conf <- def_run_4_4_nl(0.5, 0.9, lambda_h, lambda_l, alpha, model=det_non_linear, M=m_random)

non_sep <- def_run_4_4_nl(0.5, 0.5, lambda_h, lambda_l, alpha, model=det_non_linear, M=m_sep)
conf_sep <- def_run_4_4_nl(0.5, 0.9, lambda_h, lambda_l, alpha, model=det_non_linear, M=m_sep)

strat_est <- function(res, frac_h=0.5){
    return(data.frame(t=res$rr$t, rr=res$rr$CRR_h*frac_h + res$rr$CRR_l*(1-frac_h)))
}


plot_df <- rbind(non$rr %>% select(t, rr=CRR_t) %>% mutate(conf="No", mixing="Random mixing"), strat_est(conf, frac_h=0.5) %>% mutate(conf="Controlled", mixing="Random mixing"),
                conf$rr %>% select(t, rr=CRR_t) %>% mutate(conf="Yes", mixing="Random mixing"),
                    non_sep$rr %>% select(t, rr=CRR_t) %>% mutate(conf="No", mixing="Assortative mixing"), strat_est(conf_sep, frac_h=0.5) %>% mutate(conf="Controlled", mixing="Assortative mixing"),
                    conf_sep$rr %>% select(t, rr=CRR_t) %>% mutate(conf="Yes", mixing="Assortative mixing")
                    )

ggplot(plot_df) + geom_line(aes(x=t, y=rr, col=conf), size=1.7) + theme_minimal() + labs(y="Cumulative incidence ratio", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("Confounding", palette="Dark2")  + facet_wrap(.~mixing, scales="free_y")

ggsave("output/confounding.png", width=10, height=6)

m_non <-  glm(I ~ vac + offset(log(N)), data=data.frame(vac=c(0,0,1,1), risk=c(1,0,1,0), N=non$N,I=non$N - as.numeric(non$res[100, 2:5])), family="poisson")
exp(m_non$coefficients[2:3])









init.values = c(
 S = 10^3 # susceptible humans
 ) # recovered (and immune) humans
transitions = list(
  c(S = -1) # recovery
  )
init.values.SIR = c(
 S = 10^3, # susceptible humans
 I=1,
 R=0
 ) # recovered (and immune) humans
transitions.SIR = list(
  c(S = -1, I=1),
  c(I = -1, R=1) 
  )

LinearRateF <- function(x, p, t) {
 return(c(x["S"] * p$rate)) # infection rate
  
  }
NonLinearRateF <- function(x, p, t) {
 return(c(p$beta*x["S"]*x["I"]/p$N, x["I"])
        ) # infection rate

  
  }
params = list(rate=1e-2, N=1e3+1)
params_nl = list(rate=1e-1, N=1e3+1, beta=2)
r_lin <- list()
r_non <- list()

for(i in 1:1000){
    r=ssa.adaptivetau(init.values, transitions, LinearRateF, params, tf=200, tl.params=list(epsilon=0.001))
    r_lin[[i]] <- as.data.frame(r) %>% mutate(sim=i)
    r=ssa.adaptivetau(init.values.SIR, transitions.SIR, NonLinearRateF, params_nl, tf=50, tl.params=list(epsilon=0.001))
    r_non[[i]] <- as.data.frame(r) %>% mutate(sim=i)
}

res_lin <- rbindlist(r_lin)
res_non <- rbindlist(r_non)
res_lin[, cum:=10^3 - S]
res_non[, cum:=10^3 - S]
comb <- rbind(res_lin %>% mutate(type="Linear") %>% select(time, cum, type,,sim), res_non %>% mutate(type="Non-Linear")%>%select(time, cum, type, sim))
comb[, int_time:=ceiling(time)]

ggplot(comb) + geom_line(aes(x=time, y=cum, group=sim)) + facet_wrap(.~type, scales="free_x") + theme_minimal()  + scale_color_brewer("Groups", palette="Dark2") + xlab("Time") + ylab("Cumulative incidence")  + theme(text = element_text(size=20)) #+ geom_line(data=comb %>% group_by(type, int_time) %>% summarise(cum=mean(cum)), aes(x=int_time, y=cum), col="red", size=1.5)
ggsave("article/figures/stoch_ci.png", width=10, height=6, bg="white")




# Fit

N <- 100*10
pos_vac <- 10*10
pos_cont <- 20*10
t <- 10*2



NonLinearRateF <- function(x, p, t) {
    return(c(p$beta*x["S1"]*(x["I1"] + x["I2"])/p$N, x["I1"], p$alpha*p$beta*x["S2"]*(x["I1"] + x["I2"])/p$N, x["I2"])
            ) # infection rate
}

NonLinearRateFEm <- function(x, p, t) {
    return(c(p$beta*x["S1"]*(x["I1"] + x["I2"] + x["I3"])/p$N, x["I1"], p$alpha*p$beta*x["S2"]*(x["I1"] + x["I2"] + x["I3"])/p$N, x["I2"], p$beta*x["S3"]*(x["I1"] + x["I2"] + x["I3"])/p$N, x["I3"])
            ) # infection rate
}


LinearRateF <- function(x, p, t) {
    return(c(p$beta*x["S1"], 0, p$alpha*p$beta*x["S2"], 0)
            ) # infection rate
}






calc_EATEs <- function(N_vac,N_cont, pos_vac, pos_cont,I0_vac=0, I0_cont=1, t=20, n_sims=1000, inner_n=10, cores=10, epsilon=0.05, frac_cont=1, rate_func=NonLinearRateF,scale_beta=1){

    init.vals <-
        c(S1 = N_cont, # susceptible humans
            I1 = I0_cont,
            R1 = 0,
            C1 = 0,
            N1 = N_cont, 
            S2 = N_vac, # susceptible humans
            I2 = I0_vac,
            R2 = 0,
            C2 = 0,
            N2 = N_vac
        ) # recovered (and immune) humans
    transitions.vacSIR = list(
    c(S1 = -1, I1=1, C1=1),
    c(I1 = -1, R1=1),
    c(S2 = -1, I2=1, C2=1),
    c(I2 = -1, R2=1) 
    )


    if(frac_cont < 1){
        init.vals <-
        c(S1 = N_cont*frac_cont, # susceptible humans
            I1 = I0_cont*frac_cont,
            R1 = 0,
            C1 = 0,
            N1 = N_cont*frac_cont, 
            S2 = N_vac, # susceptible humans
            I2 = I0_vac,
            R2 = 0,
            C2 = 0,
            N2 = N_vac,
            S3 = N_cont*(1-frac_cont), # susceptible humans
            I3 = I0_cont*(1-frac_cont),
            R3 = 0,
            C3 = 0,
            N3 = N_cont*(1-frac_cont)

        ) # recovered (and immune) humans
    transitions.vacSIR = list(
      c(S1 = -1, I1=1, C1=1),
        c(I1 = -1, R1=1),
        c(S2 = -1, I2=1, C2=1),
     c(I2 = -1, R2=1) ,
        c(S3 = -1, I3=1, C3=1),
     c(I3 = -1, R3=1) 
        )
    }


    params_nl = list(N=N_vac + N_cont+I0_vac + I0_cont, beta=2, alpha=0.1)

    bets <- c()
    alphas <- c()
    I0s <- c()

    check_param <- function(i){
        beta <- runif(1, 1, 2.5)*scale_beta
        I0 <- sample(1:(0.1*N_cont),1)
        params_nl$beta <- beta
        alpha <- runif(1, 0.3, 0.8)
        params_nl$alpha <- alpha
        init.vals["I2"] <- I0
        if(frac_cont<1){
            init.vals["I2"] <- round(I0*frac_cont)
            init.vals["I3"] <- round(I0*(1-frac_cont))

        }
       
        r=ssa.adaptivetau(init.vals, transitions.vacSIR,rate_func, params_nl, tf=t, tl.params=list(epsilon=0.01))
        
        
        #print(N_vac - r[nrow(r), "S2"])
        #print()
       # browser()
      #  print(glue::glue("Beta: {beta}, Alpha: {alpha}, C1: {r[nrow(r), 'C1']}, C2: {r[nrow(r), 'C2']}"))
        if(abs((r[nrow(r), "C1"] - pos_cont))/(N_cont*frac_cont) < epsilon & abs((r[nrow(r), "C2"] - pos_vac))/(N_vac) < epsilon){
            return(c(beta, alpha, I0))
        } else {
            return(NULL)
        }
    }
    parallel::mclapply(1:n_sims, check_param, mc.cores=cores) -> params_found
    params_found <- params_found[!sapply(params_found, is.null)]
    betas <- sapply(params_found, function(x) x[1])
    alphas <- sapply(params_found, function(x) x[2])
    I0s <- sapply(params_found, function(x) x[3])
    
    if(length(betas) < 1){
        print("Not enough parameter sets found, consider increasing n_sims or epsilon")
        print(glue::glue("Only {length(betas)} sets found"))
        return(NULL)
    }
    #print(betas)
    print(glue::glue("Found parameter sets: {length(betas)}"))
    param_list <- lapply(1:length(betas), function(i) list(beta=betas[i], alpha=alphas[i], N=N_vac + N_cont + I0_vac + I0_cont))
    

    init.vals2 <- copy(init.vals)
    init.vals2["S1"] <- init.vals2["S1"] - 1
    init.vals2["N1"] <- init.vals2["N1"] - 1
    init.vals2["S2"] <- init.vals2["S2"] + 1
    init.vals2["N2"] <- init.vals2["N2"] + 1
        
    if(frac_cont < 1){
         init.vals <- lapply(1:length(I0s), function(i){
            tmp <- copy(init.vals)
            tmp["I2"] <- round(I0s[i]*frac_cont)
            tmp["I3"] <- round(I0s[i]*(1-frac_cont))

            return(tmp)
        })
           init.vals2 <- lapply(1:length(I0s), function(i){
            tmp <- copy(init.vals2)
            tmp["I2"] <- round(I0s[i]*frac_cont)
            tmp["I3"] <- round(I0s[i]*(1-frac_cont))

            return(tmp)
        })
    }else{
        init.vals <- lapply(1:length(I0s), function(i){
            tmp <- copy(init.vals)
            tmp["I2"]  <- I0s[i]
            return(tmp)
        })
           init.vals2 <- lapply(1:length(I0s), function(i){
            tmp <- copy(init.vals2)
            tmp["I2"]  <- I0s[i]
            return(tmp)
        })
    }
    
    eate_diff_p1 <- compare_cf(
        init.vals,
        init.vals2,
        transitions.vacSIR,
        
        param_list,
        param_list,
        rate_func,
        t,
        ((C2_Z +0.5)/N2_Z) / ((C1_X + 0.5)/N1_X),
        N=length(param_list)*inner_n,
        n_cores=cores,
        paramlist=TRUE
    
    )
    print("finished 1")
    init.vals2 <- lapply(init.vals, function(tmp){
        tmp["S1"] <- tmp["S1"] +1
        tmp["N1"] <- tmp["N1"] +1
        tmp["S2"] <- tmp["S2"] - 1
        tmp["N2"] <- tmp["N2"] - 1
        return(tmp)
    })
    eate_diff_m1 <- compare_cf(
        init.vals,
        init.vals2,
        transitions.vacSIR,
        
        param_list,
        param_list,
        rate_func,
        t,
        ((C2_Z +0.5)/N2_Z) / ((C1_X + 0.5)/N1_X),
        N=length(param_list)*inner_n,
        n_cores=cores,
        paramlist=TRUE
    )

    ind_p1 <- eate_diff_p1$diffs%>%filter(type=="ind" & time==t) %>% pull(res)
    ind_m1 <- eate_diff_m1$diffs%>%filter(type=="ind" & time==t) %>% pull(res)


    sta_p1 <- eate_diff_p1$diffs%>%filter(type=="stacked" & time==t) %>% pull(res)
    sta_m1 <- eate_diff_m1$diffs%>%filter(type=="stacked" & time==t) %>% pull(res)

    return(list(EATEs_stacked=(N_cont*frac_cont*sta_p1 + N_vac*sta_m1)/(N_cont*frac_cont + N_vac) ,EATEs_ind=(N_cont*frac_cont*ind_p1 + N_vac*ind_m1)/(N_cont*frac_cont + N_vac) , betas=betas, alphas=alphas, I0s=I0s))
}

compare_fixed_vals <- function(param_list, N_vac,N_cont, I0_vac=0, I0_cont=1, t=20, n_sims=1000, inner_n=10, cores=10, rate_func=NonLinearRateF,scale_beta=0){

    init.vals <-
        c(S1 = N_cont, # susceptible humans
            I1 = I0_cont,
            R1 = 0,
            C1 = 0,
            N1 = N_cont, 
            S2 = N_vac, # susceptible humans
            I2 = I0_vac,
            R2 = 0,
            C2 = 0,
            N2 = N_vac
        ) # recovered (and immune) humans
    transitions.vacSIR = list(
    c(S1 = -1, I1=1, C1=1),
    c(I1 = -1, R1=1),
    c(S2 = -1, I2=1, C2=1),
    c(I2 = -1, R2=1) 
    )
    init.vals2 <- copy(init.vals)
    init.vals2["S1"] <- init.vals2["S1"] - 1
    init.vals2["N1"] <- init.vals2["N1"] - 1
    init.vals2["S2"] <- init.vals2["S2"] + 1
    init.vals2["N2"] <- init.vals2["N2"] + 1
        
    print(param_list)
    eate_diff_p1 <- compare_cf(
        init.vals,
        init.vals2,
        transitions.vacSIR,
        
        param_list,
        param_list,
        rate_func,
        t,
        ((C2_Z +0.5)/N2_Z) / ((C1_X + 0.5)/N1_X),
        N=length(param_list)*inner_n,
        n_cores=cores,
        paramlist=FALSE
    
    )
    print("finished 1")
    init.vals2 <- copy(init.vals)
        init.vals2["S1"] <- init.vals2["S1"] +1
        init.vals2["N1"] <- init.vals2["N1"] +1
        init.vals2["S2"] <- init.vals2["S2"] - 1
        init.vals2["N2"] <- init.vals2["N2"] - 1

    eate_diff_m1 <- compare_cf(
        init.vals,
        init.vals2,
        transitions.vacSIR,
        
        param_list,
        param_list,
        rate_func,
        t,
        ((C2_Z +0.5)/N2_Z) / ((C1_X + 0.5)/N1_X),
        N=length(param_list)*inner_n,
        n_cores=cores,
        paramlist=FALSE
    )

    ind_p1 <- eate_diff_p1$diffs%>%filter(type=="ind" & time==t) %>% pull(res)
    ind_m1 <- eate_diff_m1$diffs%>%filter(type=="ind" & time==t) %>% pull(res)


    sta_p1 <- eate_diff_p1$diffs%>%filter(type=="stacked" & time==t) %>% pull(res)
    sta_m1 <- eate_diff_m1$diffs%>%filter(type=="stacked" & time==t) %>% pull(res)

    return(list(EATEs_stacked=(N_cont*sta_p1 + N_vac*sta_m1)/(N_cont + N_vac) ,EATEs_ind=(N_cont*ind_p1 + N_vac*ind_m1)/(N_cont + N_vac) ))
}



compare_vac_cov <- function(cov1, cov2,rate_func, params,T=20, N=1000, I0=1,n_sim=100, cores=10){
 init.vals1 <-
        c(S1 = round(N*(1-cov1)), # susceptible humans
            I1 = I0,
            R1 = 0,
            C1 = 0,
            N1 = N,
            S2 =round( N*cov1),
            I2 = 0,
            R2 = 0,
            C2 = 0,
            N2 = N
        ) # recovered (and immune) humans
 init.vals2 <-
        c(S1 = round(N*(1-cov2)), # susceptible humans
            I1 = I0,
            R1 = 0,
            C1 = 0,
            N1 = N,
            S2 = round(N*cov2),
            I2 = 0,
            R2 = 0,
            C2 = 0,
            N2 = N
        ) # recovered (and immune) humans
    transitions.vacSIR = list(
    c(S1 = -1, I1=1, C1=1),
    c(I1 = -1, R1=1),
    c(S2 = -1, I2=1, C2=1),
    c(I2 = -1, R2=1) 
    )
    
    ind1 <- parallel::mclapply(1:n_sim, function(i) regularise(as.data.frame(ssa.adaptivetau(init.vals1, transitions.vacSIR, rate_func, params, tf=T)), seq(0, T, 0.1))%>%mutate(sim=i), mc.cores = cores)
    ind2 <- parallel::mclapply(1:n_sim, function(i) regularise(as.data.frame(ssa.adaptivetau(init.vals2, transitions.vacSIR, rate_func, params, tf=T)), seq(0, T, 0.1))%>%mutate(sim=i), mc.cores = cores)

    tot <- rbindlist(lapply(1:n_sim, function(i) data.frame(time=ind1[[i]]$time, tot1=ind1[[i]]$C1 + ind1[[i]]$C2, tot2=ind2[[i]]$C1 + ind2[[i]]$C2, sim=i) %>% mutate(inc1=tot1 - lag(tot1), inc2=tot2 - lag(tot2))))

    return(tot)

}

rs <- compare_vac_cov(0.2, 0.5, NonLinearRateF,list(beta=3, alpha=0.1, N=1000),N=1000, T=30, cores=100, n_sim=1000 )



rs %>% filter(time==30) %>% summarise(sum(tot1 < 10), sum(tot2 < 10))

ggplot(rs %>% select(time, sim, tot1, tot2) %>% filter(time %in% c(1, 2, 3, 5)) %>% tidyr::pivot_longer(starts_with("tot"))) + geom_density(aes(x=value, fill=name),bw=10) + facet_grid(name~time)
ggsave("tot_diff.png")

ggplot(rs %>% filter(sim < 10) %>% select(time, sim, inc1, inc2) %>% tidyr::pivot_longer(starts_with("inc"))) + geom_line(aes(x=time, y=value, group=sim), alpha=0.5) + facet_wrap(.~name)
ggsave("timelines.png")


EATEs <- calc_EATEs(100, 100, 20, 40,I0_cont = 2, t=6, n_sims=20000, inner_n=2, cores=120, frac_cont =1, scale_beta=1.1 )

EATEs_linear <- calc_EATEs(100, 100, 20, 40,I0_cont = 2, t=6, n_sims=2000, inner_n=2, cores=100, frac_cont =1, scale_beta=0.05, rate_func=LinearRateF )

EATEs_fixed <- compare_fixed_vals(list(beta=1.6, alpha=0.47, N=200), 100, 100,I0_cont = 6, t=6, inner_n=100, cores=100)

EATEs_fixed_lin <- compare_fixed_vals(list(beta=0.084, alpha=0.49, N=200), 100, 100,I0_cont = 6, t=6, inner_n=100, cores=100, rate_func=LinearRateF)





comb <- rbind(data.frame(type="Non-Linear", ct="Ind", EATE=EATEs$EATEs_ind),
                data.frame(type="Non-Linear", ct="Stacked", EATE=EATEs$EATEs_stacked),
                data.frame(type="Linear", ct="Ind", EATE=EATEs_linear$EATEs_ind),
                data.frame(type="Linear", ct="Stacked", EATE=EATEs_linear$EATEs_stacked),
                data.frame(type="Non-Linear Fixed", ct="Ind", EATE=EATEs_fixed$EATEs_ind),
                data.frame(type="Non-Linear Fixed", ct="Stacked", EATE=EATEs_fixed$EATEs_stacked),
                data.frame(type="Linear Fixed", ct="Ind", EATE=EATEs_fixed_lin$EATEs_ind),
                data.frame(type="Linear Fixed", ct="Stacked", EATE=EATEs_fixed_lin$EATEs_stacked))



ggplot(comb) + geom_density(aes(x=EATE, fill=type), bw=0.04) + theme_minimal() + xlim(c(0,1.5)) + facet_grid(type~ct)
ggsave("dens.png")


ggplot(data.frame(x=EATEs$EATEs_stacked)) + geom_density(aes(x=x)) + theme_minimal() + xlim(c(0,1.5))
ggsave("dens.png")





quantile(EATEs$EATEs_stacked, probs=c(0.025, 0.5, 0.975))


EATEs <- calc_EATEs(50, 1000, 10, 20,I0_cont = 80, t=20, n_sims=2000, inner_n=1, cores=120, frac_cont= 0.05, rate_func=NonLinearRateFEm ,scale_beta=1)
quantile(EATEs$EATEs_ind, probs=c(0.025, 0.5, 0.975))
ggplot(E)


res <- list()
full_res <- list()
cores <- 120
for(i in c(2,3,4)){
    print(i)
    print("Non-Linear")
    EATEs <- calc_EATEs(50*i, 50*i, 10*i, 20*i,I0_cont = 5*i, t=20, n_sims=10000, inner_n=1, cores=cores, frac_cont=1 )
    
    ci_ea <- quantile(EATEs$EATEs_stacked, probs=c(0.025, 0.5, 0.975))
    res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,med=ci_ea[2], eate_l=ci_ea[1], eate_u=ci_ea[3], type="Full population", sto_type="Stacked")
    full_res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,EATE=EATEs$EATEs_stacked, sim=1:length(EATEs$EATEs_stacked), type="Full population", sto_type="Stacked")
    ci_ea <- quantile(EATEs$EATEs_ind, probs=c(0.025, 0.5, 0.975))
    res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,med=ci_ea[2], eate_l=ci_ea[1], eate_u=ci_ea[3], type="Full population", sto_type="Indivudal")
    full_res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,EATE=EATEs$EATEs_ind, sim=1:length(EATEs$EATEs_ind), type="Full population", sto_type="Individual")
    
    # print("Embeded")
    # EATEs <- calc_EATEs(50*i, 1000, 10*i, 20*i, I0_cont = 100, t= 20, n_sims=10000, inner_n=2, cores=cores, frac_cont =0.05*i, rate_func = NonLinearRateFEm, scale_beta=1.2)
    # ci_ea <- quantile(EATEs$EATEs_stacked, probs=c(0.025, 0.5, 0.975))
    # res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i, med=ci_ea[2], eate_l=ci_ea[1], eate_u=ci_ea[3], type="Embeded", sto_type="Stacked")
    # full_res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,EATE=EATEs$EATEs_stacked, sim=1:length(EATEs$EATEs_stacked), type="Embeded", sto_type="Stacked")
    # ci_ea <- quantile(EATEs$EATEs_ind, probs=c(0.025, 0.5, 0.975))
    # res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i, med=ci_ea[2], eate_l=ci_ea[1], eate_u=ci_ea[3], type="Embeded", sto_type="Indivudal")
    #     full_res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,EATE=EATEs$EATEs_ind, sim=1:length(EATEs$EATEs_stacked), type="Embeded", sto_type="Individual")
    # print("Linear")
    EATEs <- calc_EATEs(50*i, 50*i, 10*i, 20*i,I0_cont = 5*i, t=20, n_sims=5000, inner_n=4, cores=cores, frac_cont=1 , rate_func=LinearRateF, scale_beta=0.02)
    ci_ea <- quantile(EATEs$EATEs_stacked, probs=c(0.025, 0.5, 0.975))
    res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,med=ci_ea[2], eate_l=ci_ea[1], eate_u=ci_ea[3], type="Linear", sto_type="Stacked")
    full_res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,EATE=EATEs$EATEs_stacked, sim=1:length(EATEs$EATEs_stacked), type="Linear", sto_type="Stacked")
    ci_ea <- quantile(EATEs$EATEs_ind, probs=c(0.025, 0.5, 0.975))
    res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,med=ci_ea[2], eate_l=ci_ea[1], eate_u=ci_ea[3], type="Linear", sto_type="Indivudal")
    full_res[[length(res)+1]] <- data.frame(N=50*i, pos_vac=10*i, pos_cont=20*i,EATE=EATEs$EATEs_ind, sim=1:length(EATEs$EATEs_stacked), type="Linear", sto_type="Individual")
    
}

res <- rbindlist(res)
full_res <- rbindlist(full_res)


library(ggplot2)
pos_vac <- 20
pos_cont <- 40
ggplot(res) + geom_linerange(aes(x=N, y=med, ymin=eate_l, ymax=eate_u, col=type), position=position_dodge(width=5), size=1.2) + geom_point(aes(x=N, y=med, col=type), position=position_dodge(width=5), size=3) + theme_minimal() + labs(y="EATE", x="Sample size") + theme(text = element_text(size=20)) + scale_color_brewer("Method", palette="Dark2") + geom_hline(yintercept = pos_vac/pos_cont, col="red", linetype="dashed") + facet_wrap(.~sto_type)
ggsave("eate.png", width=10)


ggplot(full_res %>% filter(N==100)) + geom_density(aes(x=EATE, fill=type),alpha=0.3) + facet_wrap(.~sto_type) + xlim(c(0, 1.5))
ggsave("eate_100.png", width=12)

res$width <- res$eate_u - res$eate_l
ggplot(res) + geom_line(aes(x=N, y=width, col=type), size=1.2) + geom_point(aes(x=N, y=width, col=type), size=3) + theme_minimal() + labs(y="Width of 95% CI", x="Sample size") + theme(text = element_text(size=20)) + scale_color_brewer("Method", palette="Dark2") #+ geom_hline(yintercept = pos_vac/pos_cont, col="red", linetype="dashed") + ylim(0,5)
w
ggplot(data.frame(EATEs=EATEs, x=1:length(EATEs))) + geom_point(aes(x=x, y=EATEs)) + theme_minimal() + labs(y="EATE", x="Sample") + theme(text = element_text(size=20)) + geom_hline(yintercept = mean(EATEs), col="red", linetype="dashed")


bootstrap_ci <- function(res, n=3000, conf=0.8){
    #print(res)
    data <- data.frame(I=c(c(rep(1, res$T), rep(0, res$N_vac - res$T)), c(rep(1, res$C), rep(0, res$N_cont - res$C))),
                       group=(c(rep("T", res$N_vac), rep("C", res$N_cont))))

    if(res$C == 0 | res$T == 0){
        data[res$N/2 -1, "I"] <- 0.5
        data[res$N, "I"] <- 0.5
    }
    diff_prop <- function(data, indices) {
        d <- data[indices, ]
         treat_prop <- mean(d$I[d$group == "T"])
        control_prop <- mean(d$I[d$group == "C"])
        #if(control_prop == 0) {
            #control_prop <-  control_prop + 0.5
            #treat_prop <- treat_prop + 0.5
        #}
    return(treat_prop / control_prop)
    }

    boot_results <- boot::boot(data, diff_prop, R = n)
    ci <- boot::boot.ci(boot_results, type = "bca", conf=conf)
    

    return(rbindlist(list(data.frame(rel=boot_results$t0, rel_l=ci$bca[4], rel_u=ci$bca[5]))))
}






ci_EATEs <- quantile(EATEs, probs=c(0.025, 0.5, 0.975))

ggplot(data.frame(EATEs=EATEs) %>% filter(EATEs<10)) + geom_histogram(aes(x=EATEs), bins=30) + theme_minimal() + labs(y="Frequency", x="EATE") + theme(text = element_text(size=20)) + geom_vline(xintercept = pos_vac/pos_cont, col="red") + geom_vline(xintercept = ci$rel_l, col="red", linetype="dashed") + geom_vline(xintercept = ci$rel_u, col="red", linetype="dashed") + annotate("text", x=4, y=550, label=glue::glue("EATE={round(ci_EATEs[2],2)}, 95% CI=({round(ci_EATEs[1],2)}, {round(ci_EATEs[3],2)})"), color="blue", hjust=0, size=7) + annotate("text", x=4, y=400, label=glue::glue("Observed RR={round(pos_vac/pos_cont,2)}, 95% CI=({round(ci$rel_l,2)}, {round(ci$rel_u,2)})"), color="red", hjust=0, size=7) + theme(text = element_text(size=20))
ggsave("output/EATE_distribution.png", width=10, height=6)

hist(EATEs)

