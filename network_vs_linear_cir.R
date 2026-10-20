# CIR-vs-control-AR trajectories (paper fig 6). Compare a stochastic Pareto
# network (pa = 3) to a series of linear (exposure-only) models with increasing
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
#   - N = 200000, pa = 3, mean_k = 6
#   - 20% vaccinated with alpha_vac = 0.5
#   - Control = random unvaccinated subset of the same size as the vaccinated
#   - Seed = init_I randomly chosen nodes (independent of vac)
#
# Linear (per frailty sd):
#   - N = 10000 in each of vaccinated / unvaccinated group
#   - n_frailty = 100 equal-probability bins of a mean-1 Gamma(1/CV^2)
#     frailty at each slice's conditional mean (frailty_gamma, utils.R)
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
# Output every dust step so trajectories are as smooth as possible.
timepoints <- seq_len(as.integer(t_horizon / dt)) * dt

n_t <- length(timepoints)

# ---------------------------------------------------------------------------
# Network parameters
# ---------------------------------------------------------------------------

N_net        <- 200000
pl_alpha     <- 3
mean_k       <- 6
vac_frac     <- 0.20
beta_net     <- 1.0                   # R0-equivalent scale in homog limit
init_I       <- 20
n_rep_net    <- 100                   # memory ~ 4*N*n_rep*n_t doubles

# Sparse path: use sample_pareto_adj (O(N + m)) + run_stoch_adj_sparse
# (group-aggregate output). Required for N >> 10000; safe (but heavier)
# to use at any N.
use_sparse   <- TRUE
# Random unvac subset sized to match the vac group for a like-for-like
# comparison.
control_size <- round(vac_frac * N_net)

net_seed     <- 1
alloc_seed   <- 2
ctrl_seed    <- 3
init_seed    <- 4

# ---------------------------------------------------------------------------
# Linear + frailty parameters
# ---------------------------------------------------------------------------

N_lin_per_group <- 10000
n_frailty       <- 100
# Target CV^2 of the mean-1 gamma frailty (see build_frailty). The legend
# shows these nominal values; the binned frailty's realised CV^2 is a little
# lower and is printed while running.
frailty_cv2s    <- c(0, 0.5, 1, 2, 4)
beta_lin        <- 0.05               # gives control AR ~ 0.95 at t_horizon
n_rep_lin       <- 200

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

out_dir       <- "output/network_vs_linear_cir"
cache_dir     <- "output/network_cache"
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Mean-1 gamma frailty with target CV^2 (frailty_gamma() in utils.R:
# equal-probability bins at the conditional mean of each slice, so the mean is
# exactly 1). Returns (x, p); cv2 = 0 -> homogeneous (all 1).
build_frailty <- function(cv2, n_bins) {
  if (cv2 <= 1e-8) return(list(x = rep(1, n_bins), p = rep(1 / n_bins, n_bins)))
  fr <- frailty_gamma(shape = 1 / cv2, n = n_bins)
  list(x = fr$z, p = fr$p)
}

realised_cv2 <- function(fr) frailty_cv2(list(z = fr$x, p = fr$p))

# Cache the Pareto contact matrix (expensive to rebuild at N = 10000).
# Dense variant returns the full N x N matrix; sparse variant returns the
# adjacency-list object from sample_pareto_adj (used with
# run_stoch_adj_sparse and unblocks N ~ 10^5).
get_or_build_c_ij <- function(N, pl_alpha, mean_k, seed) {
  key  <- sprintf("c_ij_N%d_pa%s_k%d_seed%d.rds",
                  N, format(pl_alpha, trim = TRUE), mean_k, seed)
  path <- file.path(cache_dir, key)
  if (file.exists(path)) {
    message(glue("Loading cached contact matrix: {path}"))
    return(readRDS(path))
  }
  message(glue("Building Pareto network (N = {N}, pa = {pl_alpha}) ..."))
  set.seed(seed)
  c_ij <- get_conact_matrix_pl(N, alpha = pl_alpha, mean_k = mean_k)
  set.seed(NULL)
  saveRDS(c_ij, path)
  message(glue("Cached to {path}"))
  c_ij
}

get_or_build_adj <- function(N, pl_alpha, mean_k, seed) {
  key  <- sprintf("adj_N%d_pa%s_k%d_seed%d.rds",
                  N, format(pl_alpha, trim = TRUE), mean_k, seed)
  path <- file.path(cache_dir, key)
  if (file.exists(path)) {
    message(glue("Loading cached adjacency list: {path}"))
    return(readRDS(path))
  }
  message(glue("Building sparse Pareto adjacency (N = {N}, pa = {pl_alpha}) ..."))
  adj <- sample_pareto_adj(N, alpha = pl_alpha, mean_k = mean_k, seed = seed)
  saveRDS(adj, path)
  message(glue("Cached to {path}"))
  adj
}

# ---------------------------------------------------------------------------
# Network run
# ---------------------------------------------------------------------------

set.seed(alloc_seed)
vac <- sample(seq_len(N_net), round(vac_frac * N_net))
set.seed(NULL)
non_vac_all <- setdiff(seq_len(N_net), vac)

set.seed(ctrl_seed)
control_subset <- sample(non_vac_all, control_size)
set.seed(NULL)

# Random seed (independent of vac allocation) — avoids the seed-reach bias
# discussed in network_foi_cv.
set.seed(init_seed)
seeded <- sample(seq_len(N_net), init_I)
set.seed(NULL)
I_ini_vec <- integer(N_net); I_ini_vec[seeded] <- 1L

susept <- rep(1, N_net); susept[vac] <- alpha_vac

S0_vac     <- sum(1L - I_ini_vec[vac])
S0_control <- sum(1L - I_ini_vec[control_subset])

if (use_sparse) {
  # --- Sparse path: sample_pareto_adj + run_stoch_adj_sparse -------------
  # Contact structure and dust output both stay O(N * mean_k) instead of
  # O(N^2). Only aggregate S_vac, S_control totals per (t, rep) come back.
  adj <- get_or_build_adj(N_net, pl_alpha, mean_k, net_seed)
  degree_net <- adj$degree
  message(glue("Running network SIR — SPARSE (n_rep = {n_rep_net}, dt = {dt}) ..."))
  raw_net <- run_stoch_adj_sparse(
    adj,
    beta             = N_net * beta_net / mean_k,
    t                = t_horizon,
    I_ini            = I_ini_vec,
    groups           = list(vac = vac, ctrl = control_subset),
    susceptibility   = susept,
    gamma            = gamma,
    dt               = dt,
    timepoints       = timepoints,
    n_sim            = n_rep_net,
    cores            = 8)
  N_group_net <- attr(raw_net, "N_group")
  setDT(raw_net)
  S_vac_tot     <- matrix(raw_net$S_vac,  nrow = n_t, ncol = n_rep_net, byrow = TRUE)
  S_control_tot <- matrix(raw_net$S_ctrl, nrow = n_t, ncol = n_rep_net, byrow = TRUE)

  # Observed CV^2 of cumulative FOI across individuals in the CONTROL
  # subset, at the final timepoint. Uses the group sum / sum-of-squares
  # aggregates produced by stoch_mod_adj_sparse.R:
  #   mean = sum / N ; var = sumsq/N - (sum/N)^2 ; CV^2 = var / mean^2.
  # This is the "network process amplified" heterogeneity — the counterpart
  # to the degree CV^2 that ignores dynamical amplification.
  final_mask <- raw_net$time == max(raw_net$time)
  s_final    <- mean(raw_net$cumFOI_sum_ctrl[final_mask])
  s2_final   <- mean(raw_net$cumFOI_sumsq_ctrl[final_mask])
  N_ctrl_lin <- N_group_net[["ctrl"]]
  mean_cf    <- s_final  / N_ctrl_lin
  var_cf     <- s2_final / N_ctrl_lin - mean_cf^2
  cv2_foi_net <- if (mean_cf > 0) var_cf / mean_cf^2 else NA_real_
} else {
  # --- Dense path: cached full c_ij + run_stoch_adj ----------------------
  c_ij <- get_or_build_c_ij(N_net, pl_alpha, mean_k, net_seed)
  degree_net <- rowSums(c_ij)
  message(glue("Running network SIR — DENSE (n_rep = {n_rep_net}, dt = {dt}) ..."))
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
  # Column-wise reduction to avoid a [n_t, n_rep, N] intermediate.
  S_vac_tot     <- matrix(0, nrow = n_t, ncol = n_rep_net)
  S_control_tot <- matrix(0, nrow = n_t, ncol = n_rep_net)
  for (k in vac)
    S_vac_tot <- S_vac_tot + .dt_col_to_t_rep_matrix(raw_net[[paste0("S", k)]], n_t, n_rep_net)
  for (k in control_subset)
    S_control_tot <- S_control_tot + .dt_col_to_t_rep_matrix(raw_net[[paste0("S", k)]], n_t, n_rep_net)
  # Dense path doesn't store cumFOI directly, but we have per-node I so
  # we can integrate it after the fact. Compute FOI_i(t) at each t via
  # (I %*% t(c_ij)) * beta/mean_k, cumulative-sum, then CV^2 across
  # control_subset at the final timepoint (average across reps).
  n_ctrl <- length(control_subset)
  cf_at_t <- matrix(0, nrow = n_rep_net, ncol = n_ctrl)
  # Build I trajectories for control_subset only (memory-light).
  # (Skipped for simplicity — set cv2 to NA on the dense path; main
  # scientific runs use the sparse path.)
  cv2_foi_net <- NA_real_
}

# CV^2 of the realised contact-matrix out-degree distribution — the
# structural (mean-field) prediction. The dynamically observed FOI CV^2
# (cv2_foi_net, computed inside the sparse branch above) is typically
# much larger because the network process concentrates exposure on hubs
# and their neighbours — the amplification the linear-frailty analogue
# doesn't capture through degree alone.
cv2_deg_net <- var(degree_net) / mean(degree_net)^2
message(glue("Network degree CV^2 = {round(cv2_deg_net, 3)} ",
             "(mean k = {round(mean(degree_net), 2)})"))
if (!is.na(cv2_foi_net))
  message(glue("Network observed cumFOI CV^2 = {round(cv2_foi_net, 3)} ",
               "(control subset, t = {max(timepoints)})"))

ar_vac_rep     <- (S0_vac     - S_vac_tot)     / S0_vac
ar_control_rep <- (S0_control - S_control_tot) / S0_control

# Per-replicate CIR — what a single study of this size would observe.
# Drop replicates where AR_C = 0 at that t (undefined / inf ratio).
cir_rep <- ar_vac_rep / ar_control_rep
cir_rep[!is.finite(cir_rep)] <- NA_real_

# Per-replicate cumulative-hazard ratio (a.k.a. overall HR):
#   Lambda_g(t) = -log(1 - AR_g(t)) = -log(S_g(t) / S_g(0))
#   HR(t)       = Lambda_vac(t) / Lambda_ctrl(t)
# For a leaky vaccine on a homogeneous mixing pool this stays at
# alpha_vac at every t; the drift under heterogeneous exposure /
# depletion is the point of comparison against CIR.
hr_rep <- -log1p(-pmin(pmax(ar_vac_rep,     0), 1 - 1e-12)) /
          -log1p(-pmin(pmax(ar_control_rep, 0), 1 - 1e-12))
hr_rep[!is.finite(hr_rep)] <- NA_real_

net_dt <- data.table(
  t          = timepoints,
  ar_control = rowMeans(ar_control_rep, na.rm = TRUE),
  ar_vac     = rowMeans(ar_vac_rep,     na.rm = TRUE),
  cir_med    = apply(cir_rep, 1L, median,   na.rm = TRUE),
  cir_lo     = apply(cir_rep, 1L, quantile, probs = 0.25, na.rm = TRUE),
  cir_hi     = apply(cir_rep, 1L, quantile, probs = 0.75, na.rm = TRUE),
  cir_rom    = rowMeans(ar_vac_rep, na.rm = TRUE) /
                rowMeans(ar_control_rep, na.rm = TRUE),   # ratio-of-means for reference
  hr_med     = apply(hr_rep, 1L, median,   na.rm = TRUE),
  hr_lo      = apply(hr_rep, 1L, quantile, probs = 0.25, na.rm = TRUE),
  hr_hi      = apply(hr_rep, 1L, quantile, probs = 0.75, na.rm = TRUE),
  model      = sprintf("Network (pa = %s)", pl_alpha),   # CV^2s are in the log
  cv2        = if (is.na(cv2_foi_net)) cv2_deg_net else cv2_foi_net)

# Free the largest intermediate
rm(raw_net); invisible(gc(verbose = FALSE))

# ---------------------------------------------------------------------------
# Linear runs (one per frailty sd)
# ---------------------------------------------------------------------------

lin_rows <- list()
for (cv2_target in frailty_cv2s) {
  fr    <- build_frailty(cv2_target, n_frailty)
  cv2v  <- realised_cv2(fr)
  message(glue("Running linear + frailty (CV^2 target = {cv2_target}, ",
               "realised = {round(cv2v, 3)}) ..."))

  # Bins: 1..n_frailty = unvac, (n_frailty+1)..2n = vac. Each bin holds
  # 2 * N_per_group * p_bin individuals, split evenly between vac/unvac.
  n_total     <- round(2 * N_lin_per_group * fr$p)
  vac_counts  <- round(0.5 * n_total)
  N_vac_lin   <- sum(vac_counts)
  N_unvac_lin <- sum(n_total - vac_counts)
  N_groups    <- c(n_total - vac_counts, vac_counts)
  susceptibility <- c(fr$x, alpha_vac * fr$x)

  raw <- run_stoch_linear_dust(
    beta = beta_lin, N = N_groups, t = t_horizon,
    susceptibility = susceptibility,
    dt = dt, timepoints = timepoints,
    n_sim = n_rep_lin, cores = 8)
  setDT(raw)

  # Sum cumulative infections across bins per group.
  raw[, vac   := rowSums(.SD), .SDcols = paste0("C", (n_frailty + 1L):(2L * n_frailty))]
  raw[, unvac := rowSums(.SD), .SDcols = paste0("C", 1L:n_frailty)]

  ar_vac_r   <- matrix(raw$vac,   nrow = n_t, ncol = n_rep_lin, byrow = TRUE) / N_vac_lin
  ar_unvac_r <- matrix(raw$unvac, nrow = n_t, ncol = n_rep_lin, byrow = TRUE) / N_unvac_lin

  # Per-replicate CIR + cumulative-hazard ratio.
  cir_r <- ar_vac_r / ar_unvac_r
  cir_r[!is.finite(cir_r)] <- NA_real_
  hr_r  <- -log1p(-pmin(pmax(ar_vac_r,   0), 1 - 1e-12)) /
           -log1p(-pmin(pmax(ar_unvac_r, 0), 1 - 1e-12))
  hr_r[!is.finite(hr_r)] <- NA_real_

  lin_rows[[length(lin_rows) + 1L]] <- data.table(
    t          = timepoints,
    ar_control = rowMeans(ar_unvac_r, na.rm = TRUE),
    ar_vac     = rowMeans(ar_vac_r,   na.rm = TRUE),
    cir_med    = apply(cir_r, 1L, median,   na.rm = TRUE),
    cir_lo     = apply(cir_r, 1L, quantile, probs = 0.25, na.rm = TRUE),
    cir_hi     = apply(cir_r, 1L, quantile, probs = 0.75, na.rm = TRUE),
    cir_rom    = rowMeans(ar_vac_r, na.rm = TRUE) /
                  rowMeans(ar_unvac_r, na.rm = TRUE),
    hr_med     = apply(hr_r, 1L, median,   na.rm = TRUE),
    hr_lo      = apply(hr_r, 1L, quantile, probs = 0.25, na.rm = TRUE),
    hr_hi      = apply(hr_r, 1L, quantile, probs = 0.75, na.rm = TRUE),
    model      = sprintf("Linear (CV² = %g)", cv2_target),
    cv2        = cv2_target)
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

# Legend text; the network as plotmath, Network (alpha_pl = exponent).
legend_labels <- c(as.list(as.character(lin_levels)),
                   bquote("Network (\u03b1"[pl] * " = " * .(pl_alpha) * ")"))

# Linear curves in Dark2, the network in black so it stands out.
n_lin <- length(lin_levels)
pal   <- c(if (n_lin <= 8L) brewer.pal(max(n_lin, 3L), "Dark2")[seq_len(n_lin)]
           else colorRampPalette(brewer.pal(8L, "Dark2"))(n_lin),
           "black")

# Drop the very-early-epidemic region where AR_control < 0.5%. That part
# is noise-dominated (both AR_vac and AR_control are single-digit case
# counts per replicate at n_rep_lin / n_rep_net = 200/100), so the CIR
# ribbon is wildly wide and not informative about the systematic drift.
plot_dt <- all_dt[ar_control >= 0.005]

# CIR and HR share the same aesthetics + x-axis; melt into a long
# `metric` column so they can be facet-panelled in one plot.
long_dt <- rbind(
  plot_dt[, .(ar_control, model, metric = "CIR",
              y_med = cir_med, y_lo = cir_lo, y_hi = cir_hi)],
  plot_dt[, .(ar_control, model, metric = "cumulative HR",
              y_med = hr_med,  y_lo = hr_lo,  y_hi = hr_hi)]
)
long_dt[, metric := factor(metric, levels = c("CIR", "cumulative HR"),
                         labels = c("Cumulative incidence ratio (CIR)",
                                    "Cumulative hazard ratio (HR)"))]

p_combined <- ggplot(long_dt, aes(x = ar_control, y = y_med,
                                   colour = model, fill = model, group = model)) +
  geom_ribbon(aes(ymin = y_lo, ymax = y_hi),
              alpha = 0.18, colour = NA) +
  geom_path(linewidth = 1) +
  geom_hline(yintercept = alpha_vac, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(name = NULL, values = pal, labels = legend_labels) +
  scale_fill_manual(name   = NULL, values = pal, labels = legend_labels) +
  theme_minimal() +
  theme(text = element_text(size = 20), legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2),
         fill   = guide_legend(nrow = 2)) +
  labs(x = "Attack rate in control group", y = NULL)

ggsave(file.path(out_dir, "cir_hr_vs_ar_control.png"), p_combined, bg = "white",
       width = 14, height = 7, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
