# CIR-vs-control-AR trajectories. Compare a stochastic Pareto network
# (pa = 3) to a series of linear (exposure-only) models with increasing
# frailty (tagged by the CV^2 the frailty distribution induces).
#
# Both models produce a per-t trajectory in the plane
#     (AR_control(t), CIR(t) = AR_vac / AR_control).
# For a leaky vaccine with per-contact protection alpha_vac, CIR(0+) =
# alpha_vac; as attack rate rises, saturation and depletion drift the
# ratio upward. Frailty and network heterogeneity change the *speed*
# of that drift: super-exposed individuals get infected first, so the
# remaining susceptible pool is enriched for low-exposure ones, and
# CIR climbs faster/slower depending on the heterogeneity mix.
#
# Network:
#   - N = 10000, pa = 3, mean_k = 6
#   - 10% vaccinated with alpha_vac = 0.5
#   - Control = random 10% subset of the unvaccinated pool
#   - Seed = init_I randomly chosen nodes (independent of vac)
#
# Linear (per frailty sd):
#   - N = 1000 in each of vaccinated / unvaccinated group
#   - n_frailty = 100 bins with sd-controlled Beta weighting; frailty
#     values normalised to population mean 1 (matches the codebase
#     convention; see project_network_beta memory).
#   - Colour label = CV^2 of the induced frailty distribution.
#
# Usage: Rscript network_vs_linear_cir.R

library(dplyr)
library(data.table)
library(ggplot2)
library(glue)
library(RColorBrewer)

source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Common parameters
# ---------------------------------------------------------------------------

alpha_vac  <- 0.5
gamma      <- 1 / 5
# Dyadic dt = 1/4 so k*dt is fp-exact for any integer k -> dust2 happy.
dt         <- 0.25
t_horizon  <- 60
# Timepoints as integer * dt (fp-safe). Output every 4 dust steps = 1 unit.
timepoints <- seq(4L, as.integer(4 * t_horizon), 4L) * dt

n_t <- length(timepoints)

# ---------------------------------------------------------------------------
# Network parameters
# ---------------------------------------------------------------------------

N_net        <- 10000
pl_alpha     <- 3
mean_k       <- 6
vac_frac     <- 0.10
beta_net     <- 1.0                   # R0-equivalent scale in homog limit
init_I       <- 50
n_rep_net    <- 100                   # memory ~ 4*N*n_rep*n_t doubles
control_frac <- 0.10                  # random 10% of pop as control subset

net_seed     <- 1
alloc_seed   <- 2
ctrl_seed    <- 3
init_seed    <- 4

# ---------------------------------------------------------------------------
# Linear + frailty parameters
# ---------------------------------------------------------------------------

N_lin_per_group <- 1000
n_frailty       <- 100
frailty_sds     <- c(0.001, 0.05, 0.15, 0.25, 0.35, 0.45)
beta_lin        <- 0.03               # gives control AR ~ 0.83 at t_horizon
n_rep_lin       <- 200

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

out_dir <- "output/network_vs_linear_cir"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

cv2_of_frailty <- function(sd, n_frailty) {
  fr <- get_frailty(sd = sd, n = n_frailty)
  raw <- exp(2.5 * fr$x)
  fn  <- raw / sum(fr$p * raw)
  sum(fr$p * (fn - 1)^2)               # var / mean^2 with mean = 1
}

# ---------------------------------------------------------------------------
# Network run
# ---------------------------------------------------------------------------

message(glue("Building Pareto network (N = {N_net}, pa = {pl_alpha}) ..."))
set.seed(net_seed)
c_ij <- get_conact_matrix_pl(N_net, alpha = pl_alpha, mean_k = mean_k)
set.seed(NULL)

set.seed(alloc_seed)
vac <- sample(seq_len(N_net), round(vac_frac * N_net))
set.seed(NULL)
non_vac_all <- setdiff(seq_len(N_net), vac)

set.seed(ctrl_seed)
control_subset <- sample(non_vac_all, round(control_frac * N_net))
set.seed(NULL)

# Random seed (independent of vac allocation) — avoids the seed-reach bias
# discussed in network_foi_cv.
set.seed(init_seed)
seeded <- sample(seq_len(N_net), init_I)
set.seed(NULL)
I_ini_vec <- integer(N_net); I_ini_vec[seeded] <- 1L

susept <- rep(1, N_net); susept[vac] <- alpha_vac

message(glue("Running network SIR (n_rep = {n_rep_net}, dt = {dt}) ..."))
raw_net <- run_stoch_adj(
  c_ij,
  beta           = N_net * beta_net / mean_k,
  t              = t_horizon,
  I_ini          = I_ini_vec,
  susceptibility = susept,
  gamma          = gamma,
  dt             = dt,
  timepoints     = timepoints,
  n_sim          = n_rep_net,
  cores          = 8)
setDT(raw_net)

# Accumulate S_vac and S_control over reps *without* materialising the
# full [n_t, n_rep, N] array — network sim has N=10000, so this saves
# ~1 GB.
S_vac_tot     <- matrix(0, nrow = n_t, ncol = n_rep_net)
S_control_tot <- matrix(0, nrow = n_t, ncol = n_rep_net)
for (k in vac)
  S_vac_tot <- S_vac_tot + .dt_col_to_t_rep_matrix(raw_net[[paste0("S", k)]], n_t, n_rep_net)
for (k in control_subset)
  S_control_tot <- S_control_tot + .dt_col_to_t_rep_matrix(raw_net[[paste0("S", k)]], n_t, n_rep_net)

S0_vac     <- sum(1L - I_ini_vec[vac])
S0_control <- sum(1L - I_ini_vec[control_subset])
ar_vac_rep     <- (S0_vac     - S_vac_tot)     / S0_vac
ar_control_rep <- (S0_control - S_control_tot) / S0_control

# Take means across replicates first, then the ratio: E[X/Y] is a
# biased estimator of E[X]/E[Y] (Jensen), and at low AR that bias
# pushes CIR upward substantially (Y small, 1/Y convex).
ar_vac_mean     <- rowMeans(ar_vac_rep,     na.rm = TRUE)
ar_control_mean <- rowMeans(ar_control_rep, na.rm = TRUE)

net_dt <- data.table(
  t          = timepoints,
  ar_control = ar_control_mean,
  ar_vac     = ar_vac_mean,
  cir        = ar_vac_mean / ar_control_mean,
  model      = sprintf("Network (pa = %s)", pl_alpha),
  cv2        = NA_real_)

# Free the largest intermediate
rm(raw_net); invisible(gc(verbose = FALSE))

# ---------------------------------------------------------------------------
# Linear runs (one per frailty sd)
# ---------------------------------------------------------------------------

lin_rows <- list()
for (sd in frailty_sds) {
  cv2v <- cv2_of_frailty(sd, n_frailty)
  message(glue("Running linear + frailty (sd = {sd}, CV^2 = {round(cv2v, 3)}) ..."))

  out <- run_stoch_frailty_linear(
    alpha = alpha_vac, sd = sd, beta = beta_lin, f = 0.5,
    N = N_lin_per_group, t = t_horizon, n_frailty = n_frailty,
    timepoints = timepoints, n_sim = n_rep_lin, cores = 8,
    method = "dust", dt = dt)
  setDT(out)

  # Reconstruct group totals to normalise. run_stoch_frailty_linear uses
  # n_total = round(2 * N * fr$p) then vac_counts = round(0.5 * n_total).
  fr <- get_frailty(sd = sd, n = n_frailty)
  n_total     <- round(2 * N_lin_per_group * fr$p)
  vac_counts  <- round(0.5 * n_total)
  N_vac_lin   <- sum(vac_counts)
  N_unvac_lin <- sum(n_total - vac_counts)

  ar_vac_r   <- matrix(out$vac,   nrow = n_t, ncol = n_rep_lin, byrow = TRUE) / N_vac_lin
  ar_unvac_r <- matrix(out$unvac, nrow = n_t, ncol = n_rep_lin, byrow = TRUE) / N_unvac_lin

  # Ratio-of-means, not mean-of-ratios (see note above network block).
  ar_vac_mean   <- rowMeans(ar_vac_r,   na.rm = TRUE)
  ar_unvac_mean <- rowMeans(ar_unvac_r, na.rm = TRUE)

  lin_rows[[length(lin_rows) + 1L]] <- data.table(
    t          = timepoints,
    ar_control = ar_unvac_mean,
    ar_vac     = ar_vac_mean,
    cir        = ar_vac_mean / ar_unvac_mean,
    model      = sprintf("Linear (CV² = %.3f)", cv2v),
    cv2        = cv2v)
}
lin_dt <- rbindlist(lin_rows)

# ---------------------------------------------------------------------------
# Combine + plot
# ---------------------------------------------------------------------------

all_dt <- rbind(net_dt, lin_dt, fill = TRUE)
fwrite(all_dt, file.path(out_dir, "results.csv"))

# Order: linear (increasing CV^2) then network
lin_levels <- lin_dt[, .SD[1], by = model][order(cv2)]$model
model_levels <- c(as.character(lin_levels), unique(as.character(net_dt$model)))
all_dt[, model := factor(as.character(model), levels = model_levels)]

n_col <- length(model_levels)
pal   <- if (n_col <= 8L) brewer.pal(max(n_col, 3L), "Dark2")[seq_len(n_col)] else colorRampPalette(brewer.pal(8L, "Dark2"))(n_col)

p <- ggplot(all_dt, aes(x = ar_control, y = cir,
                        colour = model, group = model)) +
  geom_path(size = 1) +
  geom_hline(yintercept = alpha_vac, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(name = NULL, values = pal) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Attack rate in control group",
       y = expression(CIR == AR[vac] / AR[control]))

ggsave(file.path(out_dir, "cir_vs_ar_control.png"), p,
       width = 10, height = 7, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
