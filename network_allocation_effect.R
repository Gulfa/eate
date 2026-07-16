# Effect of vaccine allocation on measured VE for a fixed Pareto network.
# For each pl_alpha (row), materialise one contact matrix, then for each
# of n_allocs vaccine allocations:
#   - Run n_rep stochastic factual replicates.
#   - Compute the "pseudo-VE" from the cumulative incidence ratio:
#         VE_CIR(t) = 1 - CIR(t),  CIR = (C_vac/N_vac) / (C_unvac/N_unvac)
#   - Compute the proper frozen-field VE for THIS specific allocation
#     (flip each individual's vac status via cumFOI):
#         num   = sum_{i in vac}     P_fac_i  +  sum_{i not in vac} P_vac_cf_i
#         denom = sum_{i not in vac} P_fac_i  +  sum_{i in vac}     P_unvac_cf_i
#         VE_EATE(t) = 1 - num/denom
#
# Plot: rows = pl_alpha, columns = metric ("1 - CIR" vs "VE (flipping)"),
# lines coloured by allocation seed (Dark2).
#
# Usage: Rscript network_allocation_effect.R

library(dplyr)
library(data.table)
library(ggplot2)
library(glue)
library(RColorBrewer)

source("det_model.R")
source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

pl_alphas    <- c(2, 5)         # network Pareto exponents (rows in the plot)
N            <- 200
mean_k       <- 6
alpha        <- 0.5             # vaccine susceptibility
beta         <- 1.5             # transmission rate (R0 = beta/gamma)
gamma        <- 1
init_I       <- 2
t_horizon    <- 8
dt           <- 0.01
timepoints   <- seq(1, t_horizon, 1)

n_allocs     <- 10
n_rep        <- 500             # dust replicates per allocation
cores        <- 4

out_dir      <- "output/network_allocation_effect"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Per-allocation runner: returns per-time CIR and frozen-field VE for one
# fixed (c_ij, vac) pair.
# ---------------------------------------------------------------------------

one_alloc_metrics <- function(c_ij, vac, beta, alpha, N, k_mean, gamma,
                              t_horizon, dt, timepoints, init_I,
                              n_rep, cores) {
  non_vac      <- setdiff(seq_len(N), vac)
  n_t          <- length(timepoints)
  susept       <- rep(1, N); susept[vac] <- alpha
  I_ini_vec    <- c(rep(1L, init_I), rep(0L, N - init_I))

  raw <- run_stoch_adj(c_ij,
                       beta           = N * beta / k_mean,
                       t              = t_horizon,
                       I_ini          = I_ini_vec,
                       susceptibility = susept,
                       gamma          = gamma,
                       dt             = dt,
                       timepoints     = timepoints,
                       n_sim          = n_rep,
                       cores          = cores)
  setDT(raw)

  # Per-individual factual case probability (mean over reps) and I trajectory
  P_factual <- matrix(0, nrow = n_t, ncol = N)
  I_mat     <- array(0, dim = c(n_t, n_rep, N))
  for (k in seq_len(N)) {
    Ck <- .dt_col_to_t_rep_matrix(raw[[paste0("C", k)]], n_t, n_rep)
    P_factual[, k] <- rowMeans(Ck)
    I_mat[,, k]    <- .dt_col_to_t_rep_matrix(raw[[paste0("I", k)]], n_t, n_rep)
  }

  # Per-replicate cumulative FOI per individual (assumes i stays susceptible)
  cum_foi_traj <- array(0, dim = c(n_t, n_rep, N))
  for (r in seq_len(n_rep)) {
    I_traj_r <- I_mat[, r, ]
    FI_r     <- (I_traj_r %*% t(c_ij)) * (beta / k_mean)
    cum_foi_traj[, r, ] <- .cum_trapz(FI_r, timepoints)
  }

  ve_cir  <- numeric(n_t)
  ve_eate <- numeric(n_t)
  for (it in seq_len(n_t)) {
    cfi        <- cum_foi_traj[it, , ]
    P_vac_cf   <- 1 - colMeans(exp(-alpha * cfi))     # if we vaccinated i
    P_unvac_cf <- 1 - colMeans(exp(-cfi))             # if we unvaccinated i
    P_fac      <- P_factual[it, ]

    # Observed CIR: (mean vac case prob) / (mean unvac case prob)
    ve_cir[it] <- 1 - (mean(P_fac[vac]) / mean(P_fac[non_vac]))

    # Frozen-field EATE, hybrid form (this specific allocation, no re-sample)
    num   <- sum(P_fac[vac])     + sum(P_vac_cf[non_vac])
    denom <- sum(P_fac[non_vac]) + sum(P_unvac_cf[vac])
    ve_eate[it] <- 1 - num / denom
  }

  data.table(t = timepoints, VE_CIR = ve_cir, VE_EATE = ve_eate)
}

# ---------------------------------------------------------------------------
# Sweep pl_alpha x allocation_seed
# ---------------------------------------------------------------------------

all_rows <- list()
for (pa in pl_alphas) {
  message(glue("Building network pl_alpha = {pa}"))
  set.seed(pa * 1000L)
  c_ij <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
  set.seed(NULL)

  for (as_seed in seq_len(n_allocs)) {
    message(glue("  allocation {as_seed}/{n_allocs}"))
    set.seed(as_seed)
    vac <- sample(seq_len(N), round(N / 2))
    set.seed(NULL)

    m <- one_alloc_metrics(c_ij = c_ij, vac = vac,
                           beta = beta, alpha = alpha,
                           N = N, k_mean = mean_k, gamma = gamma,
                           t_horizon = t_horizon, dt = dt,
                           timepoints = timepoints, init_I = init_I,
                           n_rep = n_rep, cores = cores)
    m[, pl_alpha := pa]
    m[, alloc    := as_seed]
    all_rows[[length(all_rows) + 1L]] <- m
  }
}

results <- rbindlist(all_rows)
fwrite(results, file.path(out_dir, "results.csv"))

# ---------------------------------------------------------------------------
# Plot: rows = pl_alpha, columns = metric, lines coloured by allocation
# ---------------------------------------------------------------------------

long <- melt(results,
             id.vars       = c("t", "pl_alpha", "alloc"),
             measure.vars  = c("VE_CIR", "VE_EATE"),
             variable.name = "metric",
             value.name    = "value")
long[, metric := factor(metric,
                        levels = c("VE_CIR", "VE_EATE"),
                        labels = c("1 - CIR (pseudo-VE)",
                                   "VE (frozen-field flipping)"))]
long[, pl_alpha_lbl := factor(paste0("Pareto exp. = ", pl_alpha),
                              levels = paste0("Pareto exp. = ", pl_alphas))]

n_alloc <- length(unique(long$alloc))
pal <- if (n_alloc <= 8) brewer.pal(max(n_alloc, 3L), "Dark2")[seq_len(n_alloc)]
       else colorRampPalette(brewer.pal(8L, "Dark2"))(n_alloc)

p <- ggplot(long, aes(x = t, y = value, colour = factor(alloc), group = alloc)) +
  geom_line(size = 0.9) +
  facet_grid(pl_alpha_lbl ~ metric, scales = "free_y") +
  scale_colour_manual(name = "allocation", values = pal) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "t", y = NULL)

ggsave(file.path(out_dir, "allocation_effect.png"), p,
       width = 12, height = 4 * length(pl_alphas), dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
