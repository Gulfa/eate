library(dplyr)
library(data.table)
library(ggplot2)
library(adaptivetau)
library(ggtext)
library(latex2exp)
source("det_model.R")
source("stoch_model.R")


alpha <- 0.5
frailties <- list()
for(lambda in c(0.05, 0.1)){
    for(sd in c(0.01, 0.1, 1/sqrt(12), 0.35, 0.4)){
        frailties[[length(frailties) + 1]] <- run_frailty(lambda, alpha, sd)$sum %>% mutate(sd=round(sd,2), alpha=alpha, lambda=lambda)
    }
}
frailties <- rbindlist(frailties) 

frailties$lambda <- factor(frailties$lambda, levels=c(0.05, 0.1), labels=c("&lambda;=0.05", "&lambda;=0.1"))

q3 <- ggplot(frailties) + geom_line(aes(x=unvac/1e3, y=vac/1e3, group = paste(lambda, sd), col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Cumulative incidence vaccinated", x="Cumulative incidence control") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")
q2 <- ggplot(frailties %>% filter(t< 100)) + geom_line(aes(x=t, y=HRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Hazard rate ratio", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + facet_grid(.~lambda) + theme(strip.text = element_markdown())
q1 <- ggplot(frailties %>% filter(t< 100)) + geom_line(aes(x=t, y=CRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Cumulative incidence rate ratio", x="Time") + theme(text = element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + facet_grid(.~lambda) +  theme(strip.text = element_markdown())

cowplot::plot_grid(q1, q2,q3, ncol=1, rel_heights = c(1,1,1))

ggsave("output/attack_rate_hetero.png", width=12, height=16)




# Compare EATE and CRR for different population sizes and vaccination coverage










# Fig 2 - EATE and attack rate for different vaccination coverage


df <- rbind(run_nl(0.25) %>% mutate(frac_vac=glue::glue("25% vaccinated")),
            run_nl(0.01) %>% mutate(frac_vac=glue::glue("1% vaccinated")),
            run_nl(0.5)  %>% mutate(frac_vac=glue::glue("50% vaccinated")))

df_eate <- rbind(
  get_EATE(0.01, 3) %>% mutate(frac_vac=glue::glue("1% vaccinated")),
  get_EATE(0.25, 3) %>% mutate(frac_vac=glue::glue("25% vaccinated")),
  get_EATE(0.5,  3) %>% mutate(frac_vac=glue::glue("50% vaccinated"))
)

p1 <- ggplot(df) + geom_line(aes(x=t, y=a, col=group), size=1.7) + theme_minimal() + labs(y="Attack rate", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2") + facet_wrap(~frac_vac)

p2 <- ggplot(df_eate %>% filter(t<=50)) + geom_line(aes(x=t, y=1 - CRR, col=frac_vac), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Vaccination Coverage", palette="Set1")

cowplot::plot_grid(p1, p2, ncol=1, labels=c("A", "B"), rel_heights=c(1,0.5))

ggsave("output/attack_rate_nl.png", width=10, height=6)






# Fig 3 - Difference between CRR and EATE

dfs <- data.frame()
for(N in c(10, 100, 1000)){
    for(cov in c(0.5)){
        dfs <- rbind(dfs, get_EATE(0.5, 3, gamma=1, N=N, t=30, susceptibility = c(1,0.5)) %>% mutate(N=paste("N =", N), cov=paste("Coverage =", cov)))
    }
}
    

ggplot(dfs %>% filter(t < 50)) + geom_line(aes(x=t, y=eate_frozen/CRR, col=N), size=1.7) + theme_minimal() + labs(y="Cumulative incidence rate ratio", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Population size", palette="Set1")  + facet_grid(.~cov) + theme(strip.text = element_markdown())


plot_df <- dfs %>% filter(t < 30 & cov == "Coverage = 0.5") %>% tidyr::pivot_longer(cols=c(full, frozen, CRR), names_to="method", values_to="eate")

plot_df$method <- factor(plot_df$method, levels=c("full", "frozen", "CRR"), labels=c(TeX("$VE_{full}$"), TeX("$VE_{Frozen}$"), TeX("$VE_{CRR}$")))



random_mixing <- ggplot(plot_df) + geom_line(aes(x=t, y=1- eate, col=N), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Population size", palette="Set1")  + facet_grid(.~method) + theme(strip.text = element_markdown()) + ggtitle("Random mixing") + theme(legend.position = "bottom")
random_mixing


q1 <- ggplot(plot_df %>% filter(t < 40)) + geom_line(aes(x=t, y=1 - eate, col=N), size=1.7) + theme_minimal() + labs(y=TeX("$VE_{full}$"), x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Population size", palette="Set1")   + theme(strip.text = element_markdown(), legend.position = "none")
q2 <- ggplot(plot_df %>% filter(t < 40)) + geom_line(aes(x=t, y=1- eate_frozen, col=N), size=1.7) + theme_minimal() + labs(y=TeX("$VE_{Frozen}$"), x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Population size", palette="Set1")   + theme(strip.text = element_markdown(), legend.position = "none")
q3 <- ggplot(plot_df %>% filter(t < 40)) + geom_line(aes(x=t, y=1 - CRR, col=N), size=1.7) + theme_minimal() + labs(y=TeX("$VE_{CRR}$"), x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Population size", palette="Set1")   + theme(strip.text = element_markdown())
random_mixing <- (q1 + q2 + q3) #+ plot_layout(guides = "collect") & theme(legend.position = "bottom")


q1 <- ggplot(dfs %>% filter(t < 100)) + geom_line(aes(x=t, y=eate_frozen/CRR, col=N), size=1.7) + theme_minimal() + labs(y="EATE/Frozen", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Population size", palette="Set1")  + facet_grid(.~cov) + theme(strip.text = element_markdown())
q2 <- ggplot(dfs %>% filter(t < 100)) + geom_line(aes(x=t, y=eate/CRR, col=N), size=1.7) + theme_minimal() + labs(y="EATE/Frozen", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Population size", palette="Set1")  + facet_grid(.~cov) + theme(strip.text = element_markdown())
(q1+q2)

#Frailty
rm(frailty_eate)
for(sd in c(0.1,0.2, 0.3, 0.4)){
    f <- get_frailty_eate(0.5, sd,R=3,  n_vac=50, method = "full", N=1000)
    f$sd <- sd
    if(!exists("frailty_eate")){
        frailty_eate <- f
    } else {
        frailty_eate <- rbind(frailty_eate, f)
    }
}
#f <- get_frailty_eate(0.5, 0.4,R=2,  n_vac=20, method = "full", N=200)

frailty_eate_sum <- frailty_eate %>% group_by(t, method, sd) %>% summarise(eate= mean(num)/mean(denom))

ggplot(frailty_eate_sum %>% filter(method=="full" & t < 20 )) + geom_line(aes(x=t, y=1-eate, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")  + theme(strip.text = element_markdown(), legend.position = "bottom") + ggtitle("Frailty - Average over simulations")

ggplot(frailty_eate %>% filter(method=="CRR" & t < 20 & sd==0.1 )) + geom_line(aes(x=t, y=1-eate, group=sim, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")  + theme(strip.text = element_markdown(), legend.position = "bottom") + ggtitle("Frailty - Individual simulations")


f$method <- factor(f$method, levels=c("full", "frozen", "CRR"), labels=c(TeX("$VE_{full}$"), TeX("$VE_{Frozen}$"), TeX("$VE_{CRR}$")))

avg <- f %>% group_by(t, method) %>% summarise(eate=mean(eate))
ggplot(avg) + geom_line(aes(x=t, y=1-eate, col=method), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Method", palette="Set1")  + theme(strip.text = element_markdown(), legend.position = "bottom") + ggtitle("Frailty - Average over simulations")
frailty <- ggplot( f %>% mutate(type="Frailty") %>% filter(t<30)) + geom_line(aes(x=t, y=1-eate, group=sim), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + theme(legend.position="none") + facet_grid(.~method) + ggtitle("Frailty") + geom_line(aes(x=t, y=1-eate), data=avg, color="red", size=1.7) + theme(strip.text = element_markdown()) 

f %>% group_by(t, method) %>% summarise(eate=mean(eate), eate_marg=mean(num)/mean(denom)) %>% ggplot() + geom_line(aes(x=t, y=1-eate), color="blue", size=1.7) + geom_line(aes(x=t, y=1-eate_marg), color="red", size=1.7) + theme_minimal() + labs(y="Effect measure", x="Time") + theme(text=element_text(size=20)) + facet_wrap(.~method) #+ scale_color_brewer("Method", palette="Set1")  + theme(strip.text = element_markdown(), legend.position = "bottom") + ggtitle("Network - Average over simulations")





# ---------------------------------------------------------------------------
# Network: mean-field simulations
# ---------------------------------------------------------------------------

c_ij <- get_conact_matrix_pl(200, 50, mean_k=6)
c
cs <- list()
for(pl_alpha in c(1, 2,5, 10)){

    c <- get_eate_network(alpha=0.5, beta=2, N=200, f=0.5, t=15, pl_alpha=pl_alpha, n_vac=20, method="full", slowdown=10)
    c$pl_alpha <- pl_alpha
    cs[[length(cs)+1]] <- c
}

cs <- rbindlist(cs)

cs_sum <-cs %>% group_by(t, method, pl_alpha) %>% summarise(eate= mean(num)/mean(denom))

ggplot(cs_sum %>% filter(method=="full" & t < 20 )) + geom_line(aes(x=t, y=1-eate, col=as.factor(pl_alpha)), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")  + theme(strip.text = element_markdown(), legend.position = "bottom") + ggtitle("Frailty - Average over simulations")

ggplot(cs %>% filter(method=="CRR" & t < 20 )) + geom_line(aes(x=t, y=1-eate, group=sim ), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")  + theme(strip.text = element_markdown(), legend.position = "bottom") + ggtitle("Frailty - Individual simulations") + facet_grid(.~pl_alpha) + theme(strip.text = element_markdown())



c1 <- get_eate_network(susceptibility = c(1, 0.5), beta=2, N=100, f=0.5, t=15, pl_alpha=1.2, n_vac=20, method="both", slowdown=10)


c2 <- get_eate_network(alpha=0.5, beta=2, N=300, f=0.5, t=15, pl_alpha=1, n_vac=10, frozen_field=TRUE)

c1 %>% group_by(t) %>% summarise(eate=mean(eate))
c2 %>% group_by(t) %>% summarise(eate=mean(eate))

c1$method <- factor(c1$method, levels=c("full", "frozen", "CRR"), labels=c(TeX("$VE_{full}$"), TeX("$VE_{Frozen}$"), TeX("$VE_{CRR}$")))
avg1 <- c1 %>% group_by(t, method) %>% summarise(eate=mean(eate))
avg1 <- avg1 %>% mutate(method_short= recode(method, "VE[full]"="VE", `VE_{Frozen}`="Frozen", "VE[CRR]"="1 - Average CIR"))

ggplot(avg1) + geom_line(aes(x=t, y=1-eate, col=method), size=1.7) + theme_minimal() + labs(y="Effect measure", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Method", palette="Set1")  + theme(strip.text = element_markdown(), legend.position = "bottom") + ggtitle("Network - Average over simulations")


network <- ggplot(c1 %>% filter(method=="VE[CRR]" )) + geom_line(aes(x=t, y=1 - eate, group=sim, color="1 - CIR"), size=1.2) + theme_minimal() + labs(y="VE = 1 - CIR", x="Time") + theme(text=element_text(size=20)) + theme(legend.position="bottom")   + geom_line(aes(x=t, y=1-eate, color=method_short), data=avg1 %>% filter(method!="VE[Frozen]"), size=3) + theme(strip.text = element_markdown()) + scale_color_manual("", values=c("1 - CIR"="#40436d", "VE"="#ec7c73", "1 - Average CIR"="#61d2b2")) 
network
ggsave("output/network.png", width=12, height=8)




(random_mixing / frailty/network)
ggsave("output/network_frailty.png", width=12, height=18)

# Stochastic simulations







li <- run_linear(N_cont=100, N_vac=100, beta=0.1, alpha=0.2, t=100, n_sim=100, cores=10)

li$CRR <- li$C2 / li$C1

ggplot(li) + geom_line(aes(x=time, y=S1, group=sim), color="blue") + geom_line(aes(x=time, y=S2, group=sim), color="red") + theme_minimal() + labs(y="Number susceptible", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")
ggplot(li) + geom_line(aes(x=time, y=CRR, group=sim), color="blue") + theme_minimal() + labs(y="CRR", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")

li %>% group_by(time) %>% summarise(CRR=exp(mean(log(CRR))), CRR_marg= mean(C2)/mean(C1)) %>% ggplot() + geom_line(aes(x=time, y=CRR), color="blue") + geom_line(aes(x=time, y=CRR_marg), color="red") + theme_minimal() + labs(y="CRR", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")

res <- list()
for(sd in c(0.1, 0.2, 0.3, 0.4)){

    li_frailty <- run_stoch_frailty_linear(alpha=0.2, sd=sd,beta=0.05, N=200 , n_frailty=30, n_sim=400, cores=10)

  #  li_frailty_coupled <- run_coupled_frailty_linear(alpha=0.2, sd=sd,beta=0.05, N=200 ,vac_counts_x=0, vac_counts_z=200, n_frailty=2, n_sim=20, cores=1)

    res[[length(res)+1]] <- li_frailty %>% group_by(time) %>% summarise(CRR=exp(mean(log(CRR))), CRR_marg= mean(vac)/mean(unvac)) %>% tidyr::pivot_longer(cols=c(CRR, CRR_marg), names_to="method", values_to="CRR") %>% mutate(sd=sd)
}
res <- rbindlist(res)
res <- res %>% mutate(method= recode(method, CRR="Outer", CRR_marg="Inner")) %>% mutate(sd=paste("SD =", sd))
ggplot(res) + geom_line(aes(x=time, y=CRR, col=method), size=1.7) + theme_minimal() + labs(y="Cumulative incidence ratio", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Expectation", palette="Dark2") + facet_wrap(~sd) 
ggsave("output/frailty_stoch.png", width=12, height=8)


li_frailty %>% group_by(time) %>% summarise(CRR=mean(CRR), CRR_marg= mean(vac)/mean(unvac)) %>% ggplot() + geom_line(aes(x=time, y=CRR), color="blue") + geom_line(aes(x=time, y=CRR_marg), color="red") + theme_minimal() + labs(y="CRR", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")

ggplot(li_frailty) + geom_line(aes(x=time, y=CRR, group=sim), color="blue") + theme_minimal() + labs(y="CRR", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")
ggplot(li_frailty) + geom_line(aes(x=time, y=unvac, group=sim), color="blue") + geom_line(aes(x=time, y=vac, group=sim), color="red") + theme_minimal() + labs(y="Number susceptible", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")




sir <- run_sir(N_cont=100, N_vac=100, beta=1.5/10, alpha=0.2, t=100, n_sim=200, cores=10, gamma=1/10, I0_cont=20, I0_vac=20)

sir$CRR <- sir$C2 / sir$C1
ggplot(sir %>% filter(time > 5)) + geom_line(aes(x=time, y=CRR, group=sim), color="blue") + theme_minimal() + labs(y="CRR", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")

d <- sir %>% group_by(time) %>% summarise(CRR=mean(CRR, na.rm=TRUE), CRR_marg= mean(C2)/mean(C1)) %>% filter(time > 1)
ggplot(d) + geom_line(aes(x=time, y=CRR_marg/CRR), color="blue") 
d %>% ggplot() + geom_line(aes(x=time, y=CRR), color="blue") + geom_line(aes(x=time, y=CRR_marg), color="red") + theme_minimal() + labs(y="CRR", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")


sir %>% filter(time==80) %>% ggplot() + geom_point(aes(x=C2, y=CRR), color="black") + theme_minimal() + labs(y="CRR", x="C1", title="CRR vs C1 at time 80") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")


sir_frailty <- run_stoch_frailty_cd(alpha=0.2,sd=0.3, R=2, t=100, n_frailty = 40, n_sim=400, cores=10, I_ini_total=5, N=200)

sir_frailty %>% group_by(time) %>% summarise(CRR=mean(CRR), CRR_marg= mean(vac)/mean(unvac)) %>% ggplot() + geom_line(aes(x=time, y=CRR), color="blue") + geom_line(aes(x=time, y=CRR_marg), color="red") + theme_minimal() + labs(y="CRR", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")

ggplot(sir_frailty) + geom_line(aes(x=time, y=vac, group=sim), color="blue") + geom_line(aes(x=time, y=unvac, group=sim), color="red") + theme_minimal() + labs(y="Cumulative incidence", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")


# Test stochasticity of comparissons



res <- list()
for(n in c(200, 1000, 5000)){
    bin  <- data.frame(C1=rbinom(10000, size=n/2, prob=110/200), C2=rbinom(10000, size=n/2, prob=66/200)) 
    l <- run_stoch_linear_dust(beta=0.1, N=c(n/2,n/2), susceptibility = c(1, 0.5), t=8, timepoints=seq(1, 8, 1), n_sim=10000, cores=10, dt=0.01)
    nl_comp <- run_stoch_cd_dust(beta=1.8, matrix(rep(1, 4), nrow=2), N=c(n/2,n/2), t=8, I_ini=c(20,20), susceptibility = c(1, 0.5),
                        gamma=1,timepoints=seq(1, 8, 1),dt=0.01,  n_sim=10000, cores=10)

    l_frailty <- run_stoch_frailty_linear(alpha=0.5, sd=0.3,beta=0.1,t=8, N=n/2 , n_frailty=5, n_sim=1000, cores=10, method="dust")
    nl_frailty <- run_stoch_frailty_cd(alpha=0.5,sd=0, sd_trans=0.35, R=1.35, gamma=1, t=8, n_frailty = 5, n_sim=10000, cores=10, I_ini_total=40, N=n/2, method="dust", dt=0.01)

    


    ratios <- rbind(l %>% filter(time==8) %>% mutate(ratio=C2/C1, type="Linear"),
                bin %>% mutate(ratio=C2/C1, type="Binomial"), 
                nl_comp %>% filter(time==8) %>% mutate(ratio=C2/C1, type="SIR - comp"), 
                l_frailty %>% filter(time==8) %>% mutate(type="Linear - frailty", C1=unvac, C2=vac, ratio=C2/C1),
                nl_frailty %>% filter(time==8) %>% mutate(type="SIR - frailty", C1=unvac, C2=vac, ratio=C2/C1), 
               # indi %>% filter(time==8) %>% mutate(ratio=C2/C1, type="Individual"), 
                fill=TRUE)
    

    res[[length(res) + 1]] <- ratios %>% group_by(type) %>% summarise(C1=mean(C1), C2=mean(C2), mean_ratio=mean(ratio),sd_ratio=sd(ratio), sd_C1=sd(C1), sd_C2=sd(C2)) %>% mutate(n=n)

}

mean(nl_comp %>% filter(time==8) %>% pull(C1))
mean(nl_frailty %>% filter(time==8) %>% pull(unvac))
sd(nl_comp %>% filter(time==8) %>% pull(C1))
sd(nl_frailty %>% filter(time==8) %>% pull(unvac))


mean(nl_comp %>% filter(time==8) %>% pull(C2))
mean(nl_frailty %>% filter(time==8) %>% pull(vac))
sd(nl_comp %>% filter(time==8) %>% pull(C2))
sd(nl_frailty %>% filter(time==8) %>% pull(vac))


cor(nl_frailty %>% filter(time==8) %>% pull(unvac), nl_frailty %>% filter(time==8) %>% pull(vac))
cor(nl_comp %>% filter(time==8) %>% pull(C1), nl_comp %>% filter(time==8) %>% pull(C2))


ggplot(rbind(nl_comp %>% mutate(type="SIR - comp"), nl_frailty %>% mutate(type="SIR - frailty", C1=unvac, C2=vac), fill=TRUE)) + geom_point(aes(x=time, y=C1, col=type), position = position_dodge(width=0.3)) + theme_minimal() + labs(y="C2", x="C1") + theme(text=element_text(size=20)) + scale_color_brewer("Model", palette="Dark2")

res <- rbindlist(res)

ggplot(res) + geom_line(aes(x=n, y=sd_C2, col=type), size=1.7) + theme_minimal() + labs(y="SD of C2/C1", x="Population size") + theme(text=element_text(size=20)) + scale_color_brewer("Model", palette="Dark2") + scale_x_log10()

bin  <- data.frame(C1=rbinom(10000, size=500, prob=110/200), C2=rbinom(10000, size=500, prob=66/200)) 
l <- run_stoch_linear_dust(beta=0.1, N=c(500,500), susceptibility = c(1, 0.5), t=8, timepoints=seq(1, 8, 1), n_sim=10000, cores=10, dt=0.01)
nl_comp <- run_stoch_cd_dust(beta=1.8, matrix(rep(1, 4), nrow=2), N=c(500,500), t=8, I_ini=c(10,10), susceptibility = c(1, 0.5),
                        gamma=1,timepoints=seq(1, 8, 1),dt=0.01,  n_sim=10000, cores=10)



mm <- matrix(sample(c(0,1), size=400*400, replace=TRUE, prob=c(0.975, 0.025)), nrow=400)

mm <- get_conact_matrix_pl(400, alpha=1.2, mean_k=6)



mean(rowSums(mm))

indi <- run_stoch_ind(mm, 2.05*400/6,  c(rep(1, 10),rep(0, 190), rep(1, 10), rep(0, 190)), susceptibility = c(rep(1, 200), rep(0.5, 200)), n_sim=1000, cores=10, dt=0.01, timepoints=seq(1, 8, 1), gamma=1)
indi$C1 <- rowSums(1 - indi[,paste("S", 11:200, sep="")])
indi$C2 <- rowSums(1 - indi[,paste("S", 211:400, sep="")])

as.matrix(indi[sim==1, paste("foi", 1:200, sep="")])

mean(indi %>% filter(time==8) %>% pull(C1))
mean(nl_comp %>% filter(time==8) %>% pull(C1))

mean(indi %>% filter(time==8) %>% pull(C2))
mean(nl_comp %>% filter(time==8) %>% pull(C2))
sd(indi %>% filter(time==8) %>% pull(C2))
sd(nl_comp %>% filter(time==8) %>% pull(C2))

nl_comp %>% filter(time==8) %>% summarise(sd(C2/C1))
indi %>% filter(time==8) %>% summarise(sd(C2/C1))



hist(indi %>% filter(time==8) %>% pull(C2), breaks=20)
hist(nl_comp %>% filter(time==8) %>% pull(C2), breaks=20)




ratios <- rbind(l %>% filter(time==8) %>% mutate(ratio=C2/C1, type="Linear"),
                bin %>% mutate(ratio=C2/C1, type="Binomial"), 
                nl_comp %>% filter(time==8) %>% mutate(ratio=C2/C1, type="SIR - comp"), 
                indi %>% filter(time==8) %>% mutate(ratio=C2/C1, type="Individual"), 
                fill=TRUE)

ratios %>% group_by(type) %>% summarise(mean_ratio=mean(ratio),sd_ratio=sd(ratio))
 
ggplot(ratios) + geom_density(aes(x=ratio, fill=type), alpha=0.5, position="identity") + theme_minimal() + labs(x="C2/C1 at time 8", y="Density") + theme(text=element_text(size=20)) + scale_fill_brewer("Model", palette="Dark2")

mod_lin <- function(beta_day=0.1, susceptibility=c(1, 0.5), t=8, N=c(500,500)){
    return(list(main=data.frame(t=1:t, vac=N[2]*(1 - exp(-beta_day*(1:t)*susceptibility[2])), unvac=N[1]*(1 - exp(-beta_day*(1:t)*susceptibility[1])))))
    
}

mod_lin(beta_day=0.1, susceptibility=c(1, 0.5), t=8, N=c(200,200))

mod_sir <- purrr::partial(run_det_cd, matrix(rep(1, 4), nrow=2), N=c(500,500), I_ini=c(10,10),t=8,
                        gamma=1)
mod_sir(beta_day=1.8, susceptibility=c(1, 0.5))

cor(indi %>% filter(time==8) %>% pull(C1), indi %>% filter(time==8) %>% pull(C2))
cor(nl_comp %>% filter(time==8) %>% pull(C1), nl_comp %>% filter(time==8) %>% pull(C2))

ggplot(rbind(indi %>% filter(time==8) %>% mutate(type="individual"), nl_comp %>% filter(time==8) %>% mutate(type="SIR - comp"), fill=TRUE)) + geom_point(aes(x=C1, y=C2/C1, color=type), alpha=0.5) + theme_minimal() + labs(x="C1", y="C2") + theme(text=element_text(size=20))

cor(l %>% filter(time==8) %>% pull(C1), l %>% filter(time==8) %>% pull(C2))

fit_lin <- fit_mod_norm(mod_lin, X_cont=204, X_vac=67, beta_ini=0.1, alpha_ini=0.5, n=10000, burn_in=2000, scale=matrix(c(0.16, 0, -0.09, 0.023), nrow=2), cor_matrix = cov(l %>% filter(time==8) %>% select(C1, C2)))
fit_nl <- fit_mod_norm(mod_sir, X_cont=204, X_vac=67, beta_ini=1.8, alpha_ini=0.5, n=10000, burn_in=2000, scale=matrix(c(0.16, 0, -0.09, 0.023), nrow=2), cor_matrix = cov(nl_comp %>% filter(time==8) %>% select(C1, C2)))
fit_nl_ind <- fit_mod_norm(mod_sir, X_cont=100, X_vac=60, beta_ini=1.8, alpha_ini=0.5, n=10000, burn_in=2000, scale=matrix(c(0.16, 0, -0.09, 0.023), nrow=2), cor_matrix = cov(indi %>% filter(time==8) %>% select(C1, C2)))

#mean(fit_lin[,1])
sd(fit_lin[,2])
sd(fit_nl[,2])
sd(fit_nl_ind[,2])



a <- rnorm(100, mean=50, sd=10)
b <- rnorm(100, 2, sd=0.5)*a# rnorm(100, mean=90, sd=10)

mean(a/b)
mean(a)/mean(b)



# Fitting

   mod <- purrr::partial(run_stoch_cd_ctmc, matrix(rep(1, 4), nrow=2), N=c(500,500), t=15, I_ini=c(10,10),
                        gamma=1,timepoints=seq(1, 8, 1),dt=0.01,  n_sim=4000, cores=10)

  mod_dust <- purrr::partial(run_stoch_cd_dust, matrix(rep(1, 4), nrow=2), N=c(200,200), t=8, I_ini=c(10,10),
                        gamma=1,timepoints=seq(1, 8, 1),dt=0.01,  n_sim=6000, cores=10)

  mod_dust_net <- purrr::partial(run_stoch_network, N=1000, pl_alpha=3, t=15,
                                                  vac_frac=0.5, gamma=1, dt=0.1,
                                                  timepoints=seq(1, 15, 1), n_sim=100, cores=10)

  mod_linear <- purrr::partial(run_stoch_linear_dust, N=c(200,200), t=8, dt=0.01,
                             timepoints=seq(1, 8, 1), n_sim=1000, cores=10)
    

    mod_test <- function(beta=1, susceptibility=c(1, 0.5)){
         p1 <- beta
         p2 <- susceptibility[2] * beta
        X1 <- rbinom(4000, size=500, prob=p1)
        X2 <- rbinom(4000, size=500, prob=p2)
        return(data.frame(C1=X1, C2=X2, time=15))
    }



    pars_lin_test <- fit_mod(mod_test, X_cont=200, X_vac=100, beta_ini=200/500, alpha_ini=0.5, n=5000, burn_in=200, scale=c(0.05, 0.05))
    quantile(pars_lin_test[,1], c(0.025, 0.5, 0.975))
    quantile(rbinom(20000, size=500, prob=200/500)/500, c(0.025, 0.5, 0.975))
    mod_net <- purrr::partial(run_stoch_network, c_ij=c_ij, N=1000, pl_alpha=1, t=6, vac=vac, c_ij=c_ij, gamma=1, dt=0.01, timepoints=seq(1, 6, 1), n_sim=1000, cores=10, I_ini=20)

    pars_lin <- fit_mod(mod_linear,  X_cont=100, X_vac=20, beta_ini=0.12, alpha_ini=0.2, n=5000, burn_in=200, scale=c(0.04, 0.03))

    pars_sir <- fit_mod(mod_dust,  X_cont=100, X_vac=20, beta_ini=2.3, alpha_ini=0.2, n=10000, burn_in=400, scale=c(0.2, 0.05))
    
    
    pars_net <- fit_mod(mod_net,  X_cont=300, X_vac=60, beta_ini=2.2, alpha_ini=0.5, n=1000, burn_in=200, scale=c(0.2, 0.05))    

    plot(ts(pars_sir[,1]))
    a <- mod_linear(beta=0.12, susceptibility=c(1, 0.2))[time==8,]
    
    mean(a$C1)
    mean(a$C2)
    
    tail(run_det_cd(matrix(1, nrow=2, ncol=2), rep(2.12,6), c(500,500),   6, c(10,10), susceptibility=c(1, 0.41), gamma=1, delta_t=0.01)$full_results,1)#["C[1]" ]


    nrow(a[abs(C1 - 200) < 5 & abs(C2 - 100) < 5, .(C1, C2, time)])/2000

    
    coda::effectiveSize(pars_l)
    
     quantile(, c(0.025, 0.5, 0.975))    
    quantile(a$C1, c(0.01, 0.5, 0.975))    

    eate_func <- purrr::partial(get_EATE_linear, vac_frac=0.5, t=9, N=400)
    eate_func_nl <- purrr::partial(get_EATE, vac_frac=0.5, t=9, N=400, gamma=1, init_frac=20/400)
    eate_func_nl_net <- purrr::partial(get_eate_network, f=0.5, t=9, N=1000,init_I=20, c_ij=c_ij, n_vac=10, mc.cores=10)


    tmp <- eate_func_nl_net(beta=0.01, susceptibility=c(1, 0.41))

    plot(pars_sir)

    tail(run_mean_field(alpha=0.5, beta=2.12, N=1000, vac_frac=0.5, t=6, pl_alpha=1)$sum)
    c_ij <- get_conact_matrix_pl(1000, 2, mean_k=6)
    vac <- sample(1:1000, 500)
    
    a <- mod_net(beta=2.12, susceptibility=c(1, 0.41))


    nw$sum[time==6,]

    tail(eate_func_nl(beta=1.9, susceptibility=c(1, 0.8)))

    tail(eate_func(beta=0.01, susceptibility=c(1, 0.41)))

    res <- calc_ve_from_pars(eate_func, pars_lin, n_cores=10)
    res_nl <- calc_ve_from_pars(eate_func_nl, pars_sir, n_cores=10)

    quantile(pars_lin[,2], c(0.025, 0.5, 0.975))
    quantile(pars_sir[,2], c(0.025, 0.5, 0.975))
    sd(pars_lin[,1])/mean(pars_lin[,1])
    sd(pars_sir[,1])/mean(pars_sir[,1])
    


    res %>% filter(t==9) %>% summarise(mean=mean(full), low=quantile(full, 0.025), up=quantile(full, 0.975), sd(C1), sd(full), sd(C2))
    res_nl %>% filter(t==9) %>% summarise(mean=mean(full), low=quantile(full, 0.025), up=quantile(full, 0.975), sd(C1), sd(full), sd(C2), sd(CRR))


    res %>% filter(t==9 & C1 < 105 & C1 > 95) %>% summarise(mean=mean(full), low=quantile(full, 0.025), up=quantile(full, 0.975), sd(C1), sd(full), sd(C2))
    res_nl %>% filter(t==9 & C1 < 105 & C1 > 95) %>% summarise(mean=mean(full), low=quantile(full, 0.025), up=quantile(full, 0.975), sd(C1), sd(full), sd(C2), sd(CRR))





    x <- rbind(res_nl %>% filter(t==9)%>% mutate(method="Non-linear"), res %>% filter(t==9) %>% mutate(method="Linear"), fill=TRUE)


    ggplot(x) + geom_density(aes(x=full, fill=method), alpha=0.5) + theme_minimal() + labs(y="Density", x="VE", title="Distribution of VE estimates at time 9") + theme(text=element_text(size=20)) + scale_color_brewer("Method", palette="Set1") + scale_fill_brewer("Method", palette="Set1") + theme(strip.text = element_markdown(), legend.position = "bottom")

    x %>% group_by(method) %>% summarise(cor(C1, full), cor(C2, full), cor(CRR, full))
<

    ggplot(x) + geom_point(aes(x=C1, y=C2), color="blue") + theme_minimal() + labs(y="VE", x="C1") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2") + facet_wrap(~method)
    ggplot(x) + geom_point(aes(x=C1, y=1- full), color="blue") + theme_minimal() + labs(y="VE", x="C1") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2") + facet_wrap(~method,scales = "free")

#    quantile(res %>% filter(t==15) %>% pull(C1), c(0.025, 0.5, 0.975))
    standard_ci(x1=200, x2=100, n1=500, n2=500)

standard_ci <- function(x1,x2,n1,n2){

    p1 <- x1/n1
    p2 <- x2/n2

    return(data.frame(mean=(p2)/p1, low=(p2)/p1 * exp(-1.96*sqrt(1/x1 - 1/n1 + 1/x2 - 1/n2)), up=(p2)/p1 * exp(1.96*sqrt(1/x1 - 1/n1 + 1/x2 - 1/n2))))
}



N <- 10

r <- 0.25
alpha <- 0.5
frail <- 2

parts <- 1:N


out <- data.frame(N_vac=integer(), N_unvac=integer(), rr=double(), frail_frac=double())
for(i in 1:10000){
    vac <- sample(parts, N/2)
    vac_frail <- sum(vac <= N/2)
    vac_norm <- N/2 - vac_frail
    unvac_frail <- N/2 - vac_frail
    unvac_norm <- N/2 - vac_norm
    N_vac <- rbinom(1, size=vac_frail, prob=r*alpha*frail) + rbinom(1, size=vac_norm, prob=r*alpha)
    N_unvac <- rbinom(1, size=unvac_frail, prob=r*frail) + rbinom(1, size=unvac_norm, prob=r)
    out <- rbind(out, data.frame(N_vac=N_vac, N_unvac=N_unvac, rr=N_vac/N_unvac, frail_frac=vac_frail/(vac_frail + vac_norm)))
}

sd(out$N_vac)
sd(rbinom(10000, size=N/2, prob=r*alpha*(1+frail)/2))# / rbinom(10000, size=N/2, prob=r*(1+frail)/2))



ggplot(out) + geom_point(aes(x=N_unvac, y=N_vac), alpha=0.5) + theme_minimal() + labs(y="Number of cases in vaccinated", x="Number of cases in unvaccinated") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")


sd(out$N_unvac)
sd(out$rr)

hist(out$frail_frac)
ggplot(out) + geom_point(aes(x=frail_frac, y=rr), alpha=0.5) + theme_minimal() + labs(y="Risk ratio", x="Fraction of frail in vaccinated") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2")

s

