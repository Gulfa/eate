# CV^2 of the observed force of infection over time in stochastic
# Pareto-network SIR, for several pl_alpha values.
#
# Three curves per facet (pl_alpha):
#   Degree distribution     — CV^2(degree) from the contact matrix.
#                             Constant in t — appears as a horizontal
#                             reference line.
#   FOI: all individuals    — CV^2 of FOI_i(t) across every node
#                             (regardless of current disease state).
#   FOI: still susceptible  — same, restricted to nodes with S_i(t) = 1
#                             at that time (the pool that would
#                             actually experience the FOI).
#
# FOI on i at time t (per replicate):
#     FOI_i(t) = beta / k_mean * sum_{j: c_ij[i,j] != 0} I_j(t)
# CV^2 is computed per replicate at each t and averaged across reps.
#
# Usage: Rscript network_foi_cv.R

library(dplyr)
library(data.table)
library(ggplot2)
library(glue)
library(RColorBrewer)
library(cowplot)

source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

pl_alphas   <- c(1.5, 3, 5)     # network Pareto exponents (one facet each)
N           <- 500
mean_k      <- 6
beta        <- 1.5              # R0 = beta/gamma in homog limit
gamma       <- 1
vac_frac    <- 0.5              # fraction vaccinated
alpha_vac   <- 0.5              # vaccine susceptibility multiplier
alloc_seed  <- 1                # for reproducible vac allocation
init_I      <- 20
t_horizon   <- 10
# dust2 checks each requested time against dt with strict fp equality:
# `time / dt` must land on an integer. Non-dyadic dt (0.05, 0.1, ...)
# make e.g. 1.2 / 0.05 = 23.999...96, not 24, and the check rejects.
# Build the timepoints as `integer_steps * dt` so the same fp
# representation of dt cancels on both sides of dust's check.
dt          <- 0.05
.step_out   <- round(0.1 / dt)                       # dust steps per output
.first_step <- round(1        / dt)
.last_step  <- round(t_horizon / dt)
timepoints  <- seq(.first_step, .last_step, .step_out) * dt
n_rep       <- 200
inner_cores <- 8

out_dir <- "output/network_foi_cv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# CV^2 helper
# ---------------------------------------------------------------------------

cv2 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  m <- mean(x)
  if (m == 0) return(NA_real_)
  var(x) / m^2
}

# ---------------------------------------------------------------------------
# Sweep pl_alpha
# ---------------------------------------------------------------------------

all_rows <- list()
for (pa in pl_alphas) {
  message(glue("pl_alpha = {pa}"))
  set.seed(pa * 1001L)
  c_ij <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
  set.seed(NULL)

  # Row-sum degree (each c_ij[i, j] is 0/1; row i lists i's neighbours).
  degree  <- rowSums(c_ij)
  cv2_deg <- cv2(degree)

  # Vaccinate a fixed fraction of the population; vaccinated individuals
  # get susceptibility = alpha_vac, others = 1.
  set.seed(alloc_seed)
  vac <- sample(seq_len(N), round(vac_frac * N))
  set.seed(NULL)
  non_vac <- setdiff(seq_len(N), vac)
  susept        <- rep(1, N)
  susept[vac]   <- alpha_vac
  N_vac         <- length(vac)
  N_unvac       <- length(non_vac)

  raw <- run_stoch_adj(
    c_ij,
    beta           = N * beta / mean_k,
    t              = t_horizon,
    I_ini          = c(rep(1L, init_I), rep(0L, N - init_I)),
    susceptibility = susept,
    gamma          = gamma,
    dt             = dt,
    timepoints     = timepoints,
    n_sim          = n_rep,
    cores          = inner_cores)
  setDT(raw)

  n_t   <- length(timepoints)
  I_arr <- array(0, dim = c(n_t, n_rep, N))
  S_arr <- array(0, dim = c(n_t, n_rep, N))
  for (k in seq_len(N)) {
    I_arr[,, k] <- .dt_col_to_t_rep_matrix(raw[[paste0("I", k)]], n_t, n_rep)
    S_arr[,, k] <- .dt_col_to_t_rep_matrix(raw[[paste0("S", k)]], n_t, n_rep)
  }

  cv2_all_rep    <- matrix(NA_real_, nrow = n_t, ncol = n_rep)
  cv2_susc_rep   <- matrix(NA_real_, nrow = n_t, ncol = n_rep)
  n_susc_rep     <- matrix(0L,        nrow = n_t, ncol = n_rep)

  for (r in seq_len(n_rep)) {
    I_traj <- I_arr[, r, ]                                     # [n_t, N]
    S_traj <- S_arr[, r, ]
    FI     <- (I_traj %*% t(c_ij)) * (beta / mean_k)           # [n_t, N]
    for (it in seq_len(n_t)) {
      foi_t <- FI[it, ]
      cv2_all_rep[it, r] <- cv2(foi_t)
      sus_mask <- S_traj[it, ] > 0
      n_susc_rep[it, r] <- sum(sus_mask)
      if (n_susc_rep[it, r] >= 2L)
        cv2_susc_rep[it, r] <- cv2(foi_t[sus_mask])
    }
  }

  # Attack rate per replicate at each t = fraction of the population
  # infected. Per-group AR restricts to (un)vaccinated indices; CIR is
  # the ratio between them (vaccinated / unvaccinated), the classical
  # cumulative incidence ratio.
  S_total_all   <- apply(S_arr[, , , drop = FALSE],           c(1L, 2L), sum)
  S_total_vac   <- apply(S_arr[, , vac,     drop = FALSE],    c(1L, 2L), sum)
  S_total_unvac <- apply(S_arr[, , non_vac, drop = FALSE],    c(1L, 2L), sum)
  ar_rep        <- (N       - S_total_all)   / N              # [n_t, n_rep]
  ar_vac_rep    <- (N_vac   - S_total_vac)   / N_vac
  ar_unvac_rep  <- (N_unvac - S_total_unvac) / N_unvac
  cir_rep       <- ar_vac_rep / ar_unvac_rep                  # NaN when 0/0

  rows <- data.table(
    t              = timepoints,
    pl_alpha       = pa,
    cv2_degree     = cv2_deg,   # constant in t
    cv2_all        = rowMeans(cv2_all_rep,  na.rm = TRUE),
    cv2_unexposed  = rowMeans(cv2_susc_rep, na.rm = TRUE),
    ar             = rowMeans(ar_rep,       na.rm = TRUE),
    ar_vac         = rowMeans(ar_vac_rep,   na.rm = TRUE),
    ar_unvac       = rowMeans(ar_unvac_rep, na.rm = TRUE),
    cir            = rowMeans(cir_rep,      na.rm = TRUE),
    mean_n_susc    = rowMeans(n_susc_rep,   na.rm = TRUE))
  all_rows[[length(all_rows) + 1L]] <- rows
}

results <- rbindlist(all_rows)
fwrite(results, file.path(out_dir, "results.csv"))

# ---------------------------------------------------------------------------
# Plot: three rows stitched via cowplot
#   Row 1: AR(t)     — overall attack rate (fraction of pop infected)
#   Row 2: CIR(t)    — cumulative incidence ratio: AR_vac / AR_unvac
#   Row 3: CV^2(t)   — three lines per pl_alpha facet
# Shared column layout (pl_alpha facets), Dark2 palette for CV^2.
# ---------------------------------------------------------------------------

results[, pl_alpha_lbl := factor(sprintf("Pareto exp. = %s", format(pl_alpha, trim = TRUE)),
                                 levels = sprintf("Pareto exp. = %s",
                                                  format(pl_alphas, trim = TRUE)))]

long <- melt(results,
             id.vars       = c("t", "pl_alpha", "pl_alpha_lbl"),
             measure.vars  = c("cv2_degree", "cv2_all", "cv2_unexposed"),
             variable.name = "source",
             value.name    = "cv2")
long[, source := factor(source,
                        levels = c("cv2_degree", "cv2_all", "cv2_unexposed"),
                        labels = c("Degree distribution",
                                   "FOI: all individuals",
                                   "FOI: still susceptible"))]

pal <- brewer.pal(3L, "Dark2")

common_theme <- theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA))

p_ar <- ggplot(results, aes(x = t, y = ar)) +
  geom_line(size = 1, colour = "black") +
  facet_wrap(~ pl_alpha_lbl, scales = "free_y") +
  common_theme +
  labs(x = NULL, y = "AR(t)")

p_cir <- ggplot(results, aes(x = t, y = cir)) +
  geom_line(size = 1, colour = "black") +
  facet_wrap(~ pl_alpha_lbl, scales = "free_y") +
  common_theme +
  labs(x = NULL, y = expression(CIR(t) == AR[vac] / AR[unvac]))

p_cv2 <- ggplot(long, aes(x = t, y = cv2, colour = source, group = source)) +
  geom_line(size = 1) +
  facet_wrap(~ pl_alpha_lbl, scales = "free_y") +
  scale_colour_manual(name = NULL, values = pal) +
  common_theme +
  theme(legend.position = "bottom") +
  labs(x = "t", y = expression(CV^2))

p <- cowplot::plot_grid(p_ar, p_cir, p_cv2, ncol = 1L,
                        rel_heights = c(0.9, 0.9, 1.4),
                        align = "v", axis = "lr")

ggsave(file.path(out_dir, "cv2_foi.png"), p,
       width = 4 * length(pl_alphas) + 2, height = 12, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
