# Effect of vaccine allocation on measured VE for a fixed Pareto network.
# For each pl_alpha (row) materialise one contact matrix, pre-sample
# `n_allocs` vaccine allocations, and hand them to get_stoch_eate_network
# via its `vac_list` argument. The EATE function returns both:
#   method == "full_stoch": frozen-field hybrid EATE per allocation
#     -> VE_EATE(t) = 1 - eate
#   method == "CRR":       cumulative-incidence ratio per allocation
#     -> VE_CIR(t)  = 1 - eate  (i.e. 1 - (mean vac / mean unvac))
#
# Plot: rows = pl_alpha, columns = metric, lines coloured by allocation.

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

pl_alphas    <- c(2, 5)
N            <- 200
mean_k       <- 6
alpha        <- 0.5             # vaccine susceptibility
beta         <- 1.5             # R0 = beta / gamma in homog limit
gamma        <- 1
init_I       <- 2
t_horizon    <- 8
dt           <- 0.01
timepoints   <- seq(1, t_horizon, 1)

n_allocs     <- 10
n_rep        <- 500
inner_cores  <- 4               # dust threads per EATE call
outer_cores  <- 2               # mclapply across allocations inside EATE

out_dir      <- "output/network_allocation_effect"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Sweep pl_alpha; for each network build vac_list of n_allocs draws
# ---------------------------------------------------------------------------

all_rows <- list()
for (pa in pl_alphas) {
  message(glue("Building network pl_alpha = {pa}"))
  set.seed(pa * 1000L)
  c_ij <- get_conact_matrix_pl(N, alpha = pa, mean_k = mean_k)
  set.seed(NULL)

  vac_list <- lapply(seq_len(n_allocs), function(as_seed) {
    set.seed(as_seed); v <- sample(seq_len(N), round(N / 2)); set.seed(NULL); v
  })

  message(glue("  running get_stoch_eate_network with {n_allocs} fixed allocations"))
  res <- get_stoch_eate_network(
    beta = beta, susceptibility = c(1, alpha), f = 0.5, N = N,
    t = t_horizon, c_ij = c_ij, vac_list = vac_list,
    n_rep = n_rep, k_mean = mean_k, gamma = gamma,
    dt = dt, timepoints = timepoints, init_I = init_I,
    mc.cores = outer_cores, inner_cores = inner_cores)

  setDT(res)
  # Tag with an allocation index in the order vac_list was consumed.
  # get_stoch_eate_network preserves that order — its mclapply is over
  # seq_len(n_vac) with vac_list[[i]] — but sim ids come from runif(1).
  # Recover the ordering by taking the unique sims in first-seen order.
  sim_order <- unique(res$sim)
  res[, alloc := match(sim, sim_order)]
  res[, pl_alpha := pa]
  all_rows[[length(all_rows) + 1L]] <- res
}

results <- rbindlist(all_rows, fill = TRUE)
results[, VE := 1 - eate]
fwrite(results, file.path(out_dir, "results.csv"))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

plot_dt <- results[method %in% c("full_stoch", "CRR")]
plot_dt[, metric := factor(method,
                           levels = c("CRR", "full_stoch"),
                           labels = c("1 - CIR (pseudo-VE)",
                                      "VE (frozen-field flipping)"))]
plot_dt[, pl_alpha_lbl := factor(paste0("Pareto exp. = ", pl_alpha),
                                 levels = paste0("Pareto exp. = ", pl_alphas))]

n_alloc <- length(unique(plot_dt$alloc))
pal <- if (n_alloc <= 8) brewer.pal(max(n_alloc, 3L), "Dark2")[seq_len(n_alloc)]
       else colorRampPalette(brewer.pal(8L, "Dark2"))(n_alloc)

p <- ggplot(plot_dt, aes(x = t, y = VE,
                          colour = factor(alloc), group = alloc)) +
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
