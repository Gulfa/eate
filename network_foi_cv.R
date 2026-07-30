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
init_I      <- 5
t_horizon   <- 20
dt          <- 0.02
timepoints  <- seq(1, t_horizon, 1)
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

  raw <- run_stoch_adj(
    c_ij,
    beta           = N * beta / mean_k,
    t              = t_horizon,
    I_ini          = c(rep(1L, init_I), rep(0L, N - init_I)),
    susceptibility = rep(1, N),
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

  # Cumulative incidence rate per replicate at each t = fraction of the
  # population no longer susceptible.
  S_total  <- apply(S_arr, c(1L, 2L), sum)   # [n_t, n_rep]
  cir_rep  <- (N - S_total) / N              # [n_t, n_rep]
  cir_mean <- rowMeans(cir_rep)

  rows <- data.table(
    t              = timepoints,
    pl_alpha       = pa,
    cv2_degree     = cv2_deg,   # constant in t
    cv2_all        = rowMeans(cv2_all_rep,  na.rm = TRUE),
    cv2_unexposed  = rowMeans(cv2_susc_rep, na.rm = TRUE),
    cir            = cir_mean,
    mean_n_susc    = rowMeans(n_susc_rep,   na.rm = TRUE))
  all_rows[[length(all_rows) + 1L]] <- rows
}

results <- rbindlist(all_rows)
fwrite(results, file.path(out_dir, "results.csv"))

# ---------------------------------------------------------------------------
# Plot: two rows stitched via cowplot
#   Row 1: CIR(t)    — cumulative incidence rate, one line per pl_alpha
#                      facet (mean over replicates)
#   Row 2: CV^2(t)   — three lines per pl_alpha facet
# Shared column layout (pl_alpha facets), Dark2 palette.
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

p_cir <- ggplot(results, aes(x = t, y = cir)) +
  geom_line(size = 1, colour = "black") +
  facet_wrap(~ pl_alpha_lbl, scales = "free_y") +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = NULL, y = "CIR(t)")

p_cv2 <- ggplot(long, aes(x = t, y = cv2, colour = source, group = source)) +
  geom_line(size = 1) +
  facet_wrap(~ pl_alpha_lbl, scales = "free_y") +
  scale_colour_manual(name = NULL, values = pal) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = "t", y = expression(CV^2))

p <- cowplot::plot_grid(p_cir, p_cv2, ncol = 1L,
                        rel_heights = c(0.9, 1.4), align = "v", axis = "lr")

ggsave(file.path(out_dir, "cv2_foi.png"), p,
       width = 4 * length(pl_alphas) + 2, height = 9, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
