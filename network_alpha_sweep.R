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
alphas         <- c(0.5)                     # vaccine susceptibility values

eate_method    <- "both"                     # "full" | "frozen" | "both"
eate_slowdown  <- 1
eate_n_vac     <- 5                          # vac re-samples inside get_eate_network

# Parallelism: outer = across (pl_alpha, network_seed, alpha) cells,
# inner = across the eate_n_vac allocations inside one cell. Outer is
# the bigger axis (length(pl_alphas) * length(network_seeds) *
# length(alphas)) so default sends all cores there. Avoid nested forking
# by keeping inner_cores = 1 unless cells are few.
outer_cores    <- 4
inner_cores    <- 1

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
# Run the sweep — flat job list dispatched via mclapply
# ---------------------------------------------------------------------------

jobs <- list()
for (pl in pl_alphas) {
    for (s in network_seeds) {
        for (a in alphas) {
            jobs[[length(jobs) + 1]] <- list(pl_alpha = pl, network_seed = s, alpha = a)
        }
    }
}
n_total <- length(jobs)
message(glue("Dispatching {n_total} cells across outer_cores={outer_cores}, inner_cores={inner_cores}"))

run_cell <- function(job) {
    c_ij <- networks[[sprintf("%s|%d", format(job$pl_alpha), job$network_seed)]]
    message(glue("  cell pl_alpha={job$pl_alpha} seed={job$network_seed} alpha={job$alpha}"))
    res <- get_eate_network(beta           = beta,
                            susceptibility = c(1, job$alpha),
                            f              = vac_frac,
                            N              = N,
                            t              = t,
                            c_ij           = c_ij,
                            n_vac          = eate_n_vac,
                            method         = eate_method,
                            k_mean         = mean_k,
                            slowdown       = eate_slowdown,
                            mc.cores       = inner_cores,
                            init_I         = init_I)
    setDT(res)
    res[, pl_alpha     := job$pl_alpha]
    res[, network_seed := job$network_seed]
    res[, alpha        := job$alpha]
    res
}

results <- parallel::mclapply(jobs, run_cell,
                              mc.cores       = outer_cores,
                              mc.preschedule = FALSE)

# Surface mclapply errors instead of letting rbindlist choke on them later.
errs <- vapply(results, inherits, logical(1), what = "try-error")
if (any(errs)) {
    bad <- which(errs)
    warning(glue("{length(bad)}/{n_total} cells failed; dropping. First error:\n",
                 conditionMessage(attr(results[[bad[1]]], "condition"))))
    results <- results[!errs]
}

all_res <- rbindlist(results, fill = TRUE)
fwrite(all_res, file.path(out_dir, "raw.csv"))
message(glue("Wrote {out_dir}/raw.csv ({nrow(all_res)} rows)"))

# ---------------------------------------------------------------------------
# Summarise + plot
# ---------------------------------------------------------------------------

# VE(t) trajectory for the full method only, one row per
# (t, pl_alpha, network_seed, alpha) averaged over the eate_n_vac
# inner vac re-samples.
ve_t <- all_res[method == "full",
                .(VE      = mean(1 - eate, na.rm = TRUE),
                  VE_q025 = quantile(1 - eate, 0.025, na.rm = TRUE),
                  VE_q975 = quantile(1 - eate, 0.975, na.rm = TRUE),
                  n       = .N),
                by = .(t, pl_alpha, network_seed, alpha)]
fwrite(ve_t, file.path(out_dir, "summary.csv"))

alpha_lab <- paste(unique(ve_t$alpha), collapse = ", ")

# Trajectories: VE(t), one line per network realisation, rows = pl_alpha.
# If multiple alphas are used the lines for different alphas get mixed
# in each panel — fine for the single-alpha case the user is starting
# from; otherwise add `+ facet_grid(pl_alpha ~ alpha, ...)`.
p_lines <- ggplot(ve_t, aes(x = t, y = VE,
                            group = factor(network_seed),
                            colour = factor(network_seed))) +
    geom_line(alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    facet_wrap(~ pl_alpha, labeller = label_both, nrow = 1) +
    scale_colour_viridis_d(name = "network_seed") +
    theme_minimal(base_size = 13) +
    labs(x = "t",
         y = "VE = 1 - EATE (full)",
         title = glue("VE(t) per network — mean_k={mean_k}, N={N}, beta={beta}, alpha={alpha_lab}"))
ggsave(file.path(out_dir, "ve_trajectory.png"), p_lines,
       width = 4 * length(pl_alphas), height = 4, dpi = 130)

# Band: median + min/max across networks at each (t, pl_alpha).
band <- ve_t[, .(VE_median = median(VE),
                 VE_min    = min(VE),
                 VE_max    = max(VE)),
             by = .(t, pl_alpha, alpha)]
p_band <- ggplot(band, aes(x = t, y = VE_median,
                           group = factor(pl_alpha),
                           colour = factor(pl_alpha),
                           fill = factor(pl_alpha))) +
    geom_ribbon(aes(ymin = VE_min, ymax = VE_max), alpha = 0.2, colour = NA) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_colour_viridis_d(name = "pl_alpha") +
    scale_fill_viridis_d(name = "pl_alpha") +
    theme_minimal(base_size = 13) +
    labs(x = "t",
         y = "VE (full)",
         title = glue("Median VE(t) across {length(network_seeds)} network realisations per pl_alpha (band = min/max)"))
ggsave(file.path(out_dir, "ve_band.png"), p_band,
       width = 9, height = 5, dpi = 130)

message(glue("Wrote plots to {out_dir}/"))
message("Done.")
