# HR and CRR under susceptibility frailty in the LINEAR (constant-hazard,
# exposure-only) model, det_mod_ncd.R.
#
# Each person is infected at rate lambda * z_i (unvaccinated) or
# alpha * lambda * z_i (vaccinated), no transmission. The frailty z_i is
# discretised into bins (z_k, p_k); both arms carry N * p_k people in bin k.
# Selection on z makes the marginal HR and CRR drift away from alpha as the
# high-risk people are used up, at a rate set by CV^2(z), not by the sd of
# whatever distribution z is built from.
#
# Two frailty families:
#   frailty_beta_exp()  the original construction: z = exp(2.5 x),
#                       x ~ Beta(mean 0.5, sd) on [0,1], NOT normalised, so
#                       the mean hazard rises with sd (3.5 -> 5.3 x lambda).
#   frailty_gamma()     z ~ Gamma(shape, rate = shape): mean exactly 1,
#                       CV^2 = 1/shape, lambda is the population-mean hazard.
#
# Extracted from run_paper.R (FIG: Frailty / heterogeneity, Frailty strats),
# which sources this file.

library(dplyr)
library(data.table)
library(ggplot2)

source("utils.R")
if (!exists("det_model_ncd")) det_model_ncd <- odin::odin("det_mod_ncd.R")

# ---------------------------------------------------------------------------
# Frailty distributions: list(z = bin values, p = bin weights)
# ---------------------------------------------------------------------------

frailty_beta_exp <- function(sd, n = 200) {
  if (sd <= 0) return(list(z = 1, p = 1))
  f <- get_frailty(sd = sd, n = n)
  list(z = exp(2.5 * f$x), p = f$p)
}

# Equal-probability bins at the conditional mean of each quantile slice, so
# sum(p * z) == 1 exactly (E[Z 1{a<Z<b}] = P_{shape+1}(b) - P_{shape+1}(a)
# for a mean-1 gamma). Within-bin variance is lost, so the realised CV^2 is a
# little below 1/shape (frailty_cv2() gives the realised value; the closed-form
# check prints it).  The legend shows the nominal 1/shape.
frailty_gamma <- function(shape, n = 200) {
  if (!is.finite(shape)) return(list(z = 1, p = 1))
  q <- qgamma(seq(0, 1, length.out = n + 1), shape = shape, rate = shape)
  z <- diff(pgamma(q, shape = shape + 1, rate = shape)) * n
  list(z = z / mean(z), p = rep(1 / n, n))
}

frailty_cv2 <- function(fr) {
  m <- sum(fr$p * fr$z)
  sum(fr$p * fr$z^2) / m^2 - 1
}

# ---------------------------------------------------------------------------
# Linear model with a given frailty
# ---------------------------------------------------------------------------

run_linear_frailty <- function(lambda, alpha, frailty, N = 1000, t = 5000,
                               times = 0:t) {
  n      <- length(frailty$z)
  params <- list(n = 2 * n,
                 S_ini = N * rep(frailty$p, 2),
                 susceptibility = lambda * c(frailty$z, alpha * frailty$z))
  res <- det_model_ncd$new(user = params)$run(times)

  i_u <- seq_len(n); i_v <- n + i_u
  col <- function(nm, i) res[, paste0(nm, "[", i, "]"), drop = FALSE]

  df <- data.frame(t = res[, "t"],
                   vac   = rowSums(col("C", i_v)),
                   unvac = rowSums(col("C", i_u)),
                   inc_vac   = rowSums(col("n_SR", i_v)),
                   inc_unvac = rowSums(col("n_SR", i_u)))
  # Hazards from S itself, not N - C: once an arm is nearly exhausted, N - C
  # cancels catastrophically and the ratio blows up to noise.
  df$HRR <- (df$inc_vac / rowSums(col("S", i_v))) / (df$inc_unvac / rowSums(col("S", i_u)))
  df$CRR <- df$vac / df$unvac
  df$VE  <- 1 - df$CRR

  strats <- rbindlist(lapply(i_u, function(k) {
    d <- data.frame(t = res[, "t"], frailty_group = as.character(k),
                    vac   = res[, paste0("C[", n + k, "]")] / params$S_ini[n + k],
                    unvac = res[, paste0("C[", k, "]")]     / params$S_ini[k],
                    haz_vac   = res[, paste0("hazzard_S[", n + k, "]")],
                    haz_unvac = res[, paste0("hazzard_S[", k, "]")])
    d$HRR <- d$haz_vac / d$haz_unvac
    d$CRR <- d$vac / d$unvac
    d
  }))
  list(sum = df, strats = strats, cv2 = frailty_cv2(frailty))
}

# Original interface, kept for the existing figures.
run_frailty <- function(lambda, alpha, sd, N = 1000, t = 5000, n_frailty = 200) {
  run_linear_frailty(lambda, alpha, frailty_beta_exp(sd, n_frailty), N = N, t = t)
}

# Closed form for a mean-1 gamma frailty, as a check on the binned ODE:
#   S(t) = (1 + lambda t / k)^(-k),  HR(t) = alpha (1 + lambda t/k) / (1 + alpha lambda t/k)
gamma_closed_form <- function(lambda, alpha, shape, t) {
  Lu <- lambda * t; Lv <- alpha * Lu
  if (is.finite(shape)) {
    Su <- (1 + Lu / shape)^(-shape); Sv <- (1 + Lv / shape)^(-shape)
    HR <- alpha * (1 + Lu / shape) / (1 + Lv / shape)
  } else {
    Su <- exp(-Lu); Sv <- exp(-Lv); HR <- rep(alpha, length(t))
  }
  data.frame(t = t, HRR = HR, CRR = (1 - Sv) / (1 - Su))
}

# ---------------------------------------------------------------------------
# FIG: Frailty / heterogeneity (original beta-exp frailty)
# ---------------------------------------------------------------------------

alpha <- 0.5
N <- 1000

frailties <- list()
for (lambda in c(0.05, 0.1)) {
  for (sd in c(0.01, 0.1, 1/sqrt(12), 0.35, 0.4)) {
    frailties[[length(frailties) + 1]] <- run_frailty(lambda, alpha, sd)$sum %>%
      mutate(sd = round(sd, 2), alpha = alpha, lambda = lambda)
  }
}
frailties <- rbindlist(frailties)

frailties$lambda <- factor(frailties$lambda, levels=c(0.05, 0.1), labels=c("beta == 0.05", "beta == 0.1"))

q3 <- ggplot(frailties) + geom_line(aes(x=unvac/1e3, y=vac/1e3, group=paste(lambda, sd), col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Cumulative incidence vaccinated", x="Cumulative incidence control") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2")
q2 <- ggplot(frailties %>% filter(t< 150)) + geom_line(aes(x=t, y=HRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Hazard rate ratio", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + facet_grid(.~lambda, labeller=label_parsed)
q1 <- ggplot(frailties %>% filter(t< 150)) + geom_line(aes(x=t, y=CRR, col=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Cumulative incidence rate ratio", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + facet_grid(.~lambda, labeller=label_parsed)

cowplot::plot_grid(q1, q2, q3, ncol=1, labels=c("A", "B", "C"), rel_heights=c(1,1,1))

ggsave("output/attack_rate_hetero.png", width=12, height=16, bg="white")

# ---------------------------------------------------------------------------
# Frailty strats (original beta-exp frailty)
# ---------------------------------------------------------------------------

res <- run_frailty(0.01, 0.5, 0.4, N=1000, t=250, n_frailty=20)

ggplot(res$strats) + geom_line(aes(x=unvac, y=vac, color=frailty_group))

ggplot(res$strats) + geom_line(aes(x=unvac, y=vac/unvac, color=frailty_group)) + geom_line(aes(x=unvac/1000, y=vac/unvac), data=res$sum) + theme_minimal() + labs(y="Cumulative incidence ratio", x="Cumulative incidence control") + theme(text=element_text(size=20)) + scale_color_brewer("Frailty group", palette="Dark2") + theme(legend.position="bottom")

ggplot(res$strats) + geom_line(aes(x=t, y=HRR, color=frailty_group)) + geom_line(aes(x=t, y=HRR), data=res$sum)

res2 <- run_frailty(0.01, 0.5, 0.15, N=1000, t=250, n_frailty=20)
res3 <- run_frailty(0.01, 0.5, 0.25, N=1000, t=250, n_frailty=20)
res4 <- run_frailty(0.01, 0.5, 0.0001, N=1000, t=250, n_frailty=3)
comb <- rbind(res$sum %>% mutate(sd=0.4), res2$sum %>% mutate(sd=0.15), res3$sum %>% mutate(sd=0.25), res4$sum %>% mutate(sd=0.00001))
ggplot(comb) + geom_line(aes(x=t, y=HRR, color=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Marginal Hazard rate ratio", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + theme(legend.position="bottom")
ggsave("output/hazard_rate_hetero.png", width=10, height=6, bg="white")
ggplot(comb) + geom_line(aes(x=unvac/1000, y=CRR, color=as.factor(sd)), size=1.7) + theme_minimal() + labs(y="Marginal Cumulative incidence rate ratio", x="Fraction infected in control") + theme(text=element_text(size=20)) + scale_color_brewer("SD", palette="Dark2") + theme(legend.position="bottom")
ggsave("output/cumulative_incidence_hetero.png", width=10, height=6, bg="white")

# ---------------------------------------------------------------------------
# Same figure with a mean-1 gamma frailty, varying the shape
# ---------------------------------------------------------------------------
# shape = Inf is the homogeneous model; CV^2 = 1/shape = 0, 0.25, 0.5, 1, 2.
# Because the mean is 1, lambda is the true mean hazard and every curve has
# the same initial speed -- differences are heterogeneity alone.

gamma_shapes <- c(Inf, 4, 2, 1, 0.5)
gamma_lambdas <- c(0.05, 0.1)
# Fine early grid (so panel C starts at the origin) and a log-spaced tail out
# to t_gamma: at CV^2 = 2 the control arm reaches 99% only at
# lambda t ~ 5e3, so the run has to be long.
t_gamma <- 2e5
times_gamma <- sort(unique(c(seq(0, 10, by = 0.05),
                             exp(seq(log(10), log(t_gamma), length.out = 400)))))

gam <- list(); gam_check <- list()
for (lambda in gamma_lambdas) {
  for (k in gamma_shapes) {
    r <- run_linear_frailty(lambda, alpha, frailty_gamma(k), N = N, times = times_gamma)
    lab <- format(1 / k)
    gam[[length(gam) + 1]] <- r$sum %>% mutate(shape = k, label = lab, lambda = lambda)
    # Compare only while < 99% of the control arm is infected: beyond that the
    # survivors are the lowest bin alone and the binned HR falls back towards
    # alpha, a discretisation tail rather than anything plotted.
    cf <- gamma_closed_form(lambda, alpha, k, r$sum$t)
    ok <- r$sum$t > 0 & r$sum$unvac / N < 0.99
    gam_check[[length(gam_check) + 1]] <- data.frame(
      lambda = lambda, shape = k, cv2_realised = r$cv2,
      max_dHRR = max(abs(r$sum$HRR - cf$HRR)[ok]),
      max_dCRR = max(abs(r$sum$CRR - cf$CRR)[ok]))
  }
}
gam <- rbindlist(gam)
gam$label  <- factor(gam$label, levels = unique(gam$label))
gam$lambda <- factor(gam$lambda, levels = gamma_lambdas, labels = paste("beta ==", gamma_lambdas))

message("Binned ODE vs closed form (gamma frailty):")
print(rbindlist(gam_check), digits = 3)

g3 <- ggplot(gam %>% filter(t > 0)) + geom_line(aes(x=unvac/1e3, y=CRR, group=paste(lambda, label), col=label), size=1.7) + theme_minimal() + labs(y="CIR", x="Cumulative incidence control") + theme(text=element_text(size=20)) + scale_color_brewer(expression(CV^2), palette="Dark2")
g2 <- ggplot(gam %>% filter(t < 150)) + geom_line(aes(x=t, y=HRR, col=label), size=1.7) + theme_minimal() + labs(y="Hazard rate ratio", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer(expression(CV^2), palette="Dark2") + facet_grid(.~lambda, labeller=label_parsed)
g1 <- ggplot(gam %>% filter(t < 150)) + geom_line(aes(x=t, y=CRR, col=label), size=1.7) + theme_minimal() + labs(y="Cumulative incidence rate ratio", x="Time") + theme(text=element_text(size=20)) + scale_color_brewer(expression(CV^2), palette="Dark2") + facet_grid(.~lambda, labeller=label_parsed)

cowplot::plot_grid(g1, g2, g3, ncol=1, labels=c("A", "B", "C"), rel_heights=c(1,1,1))

ggsave("output/attack_rate_hetero_gamma.png", width=12, height=16, bg="white")
