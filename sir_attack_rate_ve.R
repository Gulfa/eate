# Attack rate per arm and VE = 1 - CIR over time in the deterministic SIR
# (det_mod_cd.R, random mixing) with a fixed leaky vaccine, alpha = 0.3,
# beta = 6, gamma = 1, at 1%, 25% and 50% coverage.
#
#   A: cumulative incidence in the control and vaccinated groups
#   B: VE = 1 - CIR from the same runs; starts at 1 - alpha and falls as the
#      control group is depleted, more so the larger the epidemic
#
# Extracted from run_paper.R (Non-linear / network effects), which sources
# this file.

library(dplyr)
library(ggplot2)

if (!exists("run_nl")) source("det_model.R")

beta  <- 6
t_max <- 20

df <- rbind(run_nl(0.25, beta=beta, t=t_max, delta_t=0.1) %>% mutate(frac_vac=glue::glue("25% vaccinated")),
            run_nl(0.01, beta=beta, t=t_max, delta_t=0.1) %>% mutate(frac_vac=glue::glue("1% vaccinated")),
            run_nl(0.5,  beta=beta, t=t_max, delta_t=0.1) %>% mutate(frac_vac=glue::glue("50% vaccinated")))

# VE = 1 - CIR, from the same runs as panel A
df_ve <- df %>% tidyr::pivot_wider(names_from=group, values_from=a) %>%
  filter(t > 0) %>% mutate(VE = 1 - Vaccinated/Control)

p1 <- ggplot(df) + geom_line(aes(x=t, y=a, col=group), size=1.7) + theme_minimal() + labs(y="Attack rate", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Group", palette="Dark2") + facet_wrap(~frac_vac)

p2 <- ggplot(df_ve) + geom_line(aes(x=t, y=VE, col=frac_vac), size=1.7) + theme_minimal() + labs(y="VE", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("Vaccination Coverage", palette="Set1")

cowplot::plot_grid(p1, p2, ncol=1, labels=c("A", "B"), rel_heights=c(1,0.5))

ggsave("output/attack_rate_nl.png", width=10, height=6, bg="white")
