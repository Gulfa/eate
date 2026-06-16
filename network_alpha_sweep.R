# Deterministic VE sweep over (pl_alpha, network_seed, alpha). For each
# pl_alpha (Pareto exponent of the contact distribution) we materialise
# `network_seeds` independent contact graphs, and for every alpha
# (vaccine-susceptibility) value compute EATE / VE through
# get_eate_network directly (no MCMC fitting). Useful for seeing how VE
# changes across degree heterogeneity, network realisations, and
# vaccine effectiveness.
#
# Usage:
#   Rscript network_alpha_sweep.R

library(data.table)
library(dplyr)
library(ggplot2)
library(glue)

source("det_model.R")
source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Parameters (edit here)
# ---------------------------------------------------------------------------

pl_alphas      <- c(2.5, 3, 5)               # Pareto exponents to sweep
mean_k         <- 6                          # mean degree
N              <- 200                        # population size
t              <- 8                          # horizon (matches run_array)
beta           <- 0.5                        # transmission rate
gamma          <- 1                          # project convention; passed via run_mean_field
init_I         <- 2                          # seed infections per run
vac_frac       <- 0.5                        # coverage

network_seeds  <- 1:10                       # contact-matrix realisations per pl_alpha
alphas         <- c(0.1, 0.25, 0.5, 0.75, 0.9)  # vaccine susceptibility values

eate_method    <- "both"                     # "full" | "frozen" | "both"
eate_slowdown  <- 1
eate_n_vac     <- 5                          # vac re-samples inside get_eate_network
mc_cores       <- 4

out_dir        <- "output/network_alpha_sweep"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Materialise networks: one per (pl_alpha, network_seed). Using a
# (pl_alpha, seed) keyed RNG so seeds are independent across pl_alphas.
# ---------------------------------------------------------------------------

message(glue("Building {length(pl_alphas)} x {length(network_seeds)} contact matrices ",
             "(N={N}, mean_k={mean_k})"))
networks <- list()
for (pl in pl_alphas) {
    for (s in network_seeds) {
        set.seed(1000L * round(pl * 10) + s)   # avoid seed collisions across pl_alphas
        key <- sprintf("%s|%d", format(pl), s)
        networks[[key]] <- get_conact_matrix_pl(N, alpha = pl, mean_k = mean_k)
    }
}
set.seed(NULL)

# ---------------------------------------------------------------------------
# Run the sweep
# ---------------------------------------------------------------------------

n_total <- length(pl_alphas) * length(network_seeds) * length(alphas)
counter <- 0
results <- list()
for (pl in pl_alphas) {
    for (s in network_seeds) {
        c_ij <- networks[[sprintf("%s|%d", format(pl), s)]]
        for (a in alphas) {
            counter <- counter + 1
            message(glue("[{counter}/{n_total}] pl_alpha={pl}  network_seed={s}  alpha={a}"))
            res <- get_eate_network(beta           = beta,
                                    susceptibility = c(1, a),
                                    f              = vac_frac,
                                    N              = N,
                                    t              = t,
                                    c_ij           = c_ij,
                                    n_vac          = eate_n_vac,
                                    method         = eate_method,
                                    k_mean         = mean_k,
                                    slowdown       = eate_slowdown,
                                    mc.cores       = mc_cores,
                                    init_I         = init_I)
            setDT(res)
            res[, pl_alpha     := pl]
            res[, network_seed := s]
            res[, alpha        := a]
            results[[length(results) + 1]] <- res
        }
    }
}

all_res <- rbindlist(results, fill = TRUE)
fwrite(all_res, file.path(out_dir, "raw.csv"))
message(glue("Wrote {out_dir}/raw.csv ({nrow(all_res)} rows)"))

# ---------------------------------------------------------------------------
# Summarise + plot
# ---------------------------------------------------------------------------

# VE = 1 - eate at the final time, averaged over the eate_n_vac
# re-samples (`sim` is the inner vac-iteration id from get_eate_network).
ve <- all_res[t == max(t),
              .(VE      = mean(1 - eate, na.rm = TRUE),
                VE_q025 = quantile(1 - eate, 0.025, na.rm = TRUE),
                VE_q975 = quantile(1 - eate, 0.975, na.rm = TRUE),
                n       = .N),
              by = .(pl_alpha, network_seed, alpha, method)]
fwrite(ve, file.path(out_dir, "summary.csv"))

# Lines: VE vs alpha, one line per network realisation, rows = pl_alpha,
# columns = method.
p_lines <- ggplot(ve, aes(x = alpha, y = VE,
                          group = factor(network_seed),
                          colour = factor(network_seed))) +
    geom_line(alpha = 0.7) +
    geom_point() +
    facet_grid(pl_alpha ~ method, labeller = label_both) +
    scale_colour_viridis_d(name = "network_seed") +
    theme_minimal(base_size = 13) +
    labs(x = "alpha (vaccine susceptibility)",
         y = "VE = 1 - EATE",
         title = glue("Network VE sweep — mean_k={mean_k}, N={N}, beta={beta}"))
ggsave(file.path(out_dir, "ve_by_alpha_lines.png"), p_lines,
       width = 10, height = 3 * length(pl_alphas), dpi = 130)

# Band: median + min/max across networks at each (pl_alpha, alpha).
band <- ve[, .(VE_median = median(VE),
               VE_min    = min(VE),
               VE_max    = max(VE)),
           by = .(pl_alpha, alpha, method)]
p_band <- ggplot(band, aes(x = alpha, y = VE_median,
                           group = factor(pl_alpha),
                           colour = factor(pl_alpha),
                           fill = factor(pl_alpha))) +
    geom_ribbon(aes(ymin = VE_min, ymax = VE_max), alpha = 0.2, colour = NA) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~ method) +
    scale_colour_viridis_d(name = "pl_alpha") +
    scale_fill_viridis_d(name = "pl_alpha") +
    theme_minimal(base_size = 13) +
    labs(x = "alpha",
         y = "VE",
         title = glue("Median VE across {length(network_seeds)} network realisations per pl_alpha (band = min/max)"))
ggsave(file.path(out_dir, "ve_band.png"), p_band,
       width = 10, height = 5, dpi = 130)

message(glue("Wrote plots to {out_dir}/"))
message("Done.")
