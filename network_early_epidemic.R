# Early-epidemic exploration on Pareto contact networks. Fixed R0, no
# vaccination. Vary the Pareto exponent `pl_alpha` and for each exponent
# run several independent network realisations; measure how the early
# growth from 20 seeded infectives up to ~500 new cases depends on both
# network structure (systematic) and specific graph realisation
# (structural stochasticity) plus within-graph transmission noise.
#
# Uses the sparse adjacency sampler + aggregate-output odin model so
# memory stays trivial and we can bump N later without re-plumbing.

library(dplyr)
library(data.table)
library(ggplot2)
library(glue)
library(RColorBrewer)

source("stoch_model.R")
source("utils.R")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

N            <- 10000
mean_k       <- 6
init_I       <- 20
gamma        <- 1 / 5             # infectious period 5 (time units)
R0           <- 2
beta_net     <- R0 * gamma        # user-facing beta; R0 = beta/gamma in the
                                  # homogeneous-network limit (see project
                                  # memory on the beta convention).

pl_alphas    <- c(1.5, 5)
n_networks   <- 5                 # independent contact-matrix realisations
                                  # per pl_alpha
n_rep        <- 100               # dust replicates per (pl_alpha, network)

t_horizon    <- 30                # enough to reach ~500 cases even on the
                                  # slower networks at R0 = 2
dt           <- 0.25              # dyadic (1/4), fp-safe with 1-step output
timepoints   <- seq_len(as.integer(t_horizon / dt)) * dt
n_t          <- length(timepoints)

target_cases <- 500               # "the first few weeks until we reach ~500"

# Independent RNG streams for network build, initial-seed choice, and dust
net_seed_offset  <- 10000L
seed_seed_offset <- 20000L

out_dir   <- "output/network_early_epidemic"
cache_dir <- "output/network_cache"
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Cache the Pareto adjacency (same cache as network_vs_linear_cir.R;
# keyed by N, pl_alpha, mean_k, seed).
# ---------------------------------------------------------------------------

get_or_build_adj <- function(N, pl_alpha, mean_k, seed) {
  key  <- sprintf("adj_N%d_pa%s_k%d_seed%d.rds",
                  N, format(pl_alpha, trim = TRUE), mean_k, seed)
  path <- file.path(cache_dir, key)
  if (file.exists(path)) {
    message(glue("Loading cached adjacency: {path}"))
    return(readRDS(path))
  }
  message(glue("Building sparse Pareto adjacency (N = {N}, pa = {pl_alpha}, seed = {seed}) ..."))
  adj <- sample_pareto_adj(N, alpha = pl_alpha, mean_k = mean_k, seed = seed)
  saveRDS(adj, path)
  message(glue("Cached to {path}"))
  adj
}

# ---------------------------------------------------------------------------
# Sweep pl_alpha x network_seed
# ---------------------------------------------------------------------------

all_rows <- list()
for (pa in pl_alphas) {
  for (nseed in seq_len(n_networks)) {
    adj_seed  <- net_seed_offset  * round(pa * 10) + nseed
    init_seed <- seed_seed_offset * round(pa * 10) + nseed

    adj <- get_or_build_adj(N, pa, mean_k, seed = adj_seed)

    cv2_deg <- var(adj$degree) / mean(adj$degree)^2
    message(glue("  degree CV² = {round(cv2_deg, 3)} ",
                 "(mean k = {round(mean(adj$degree), 2)})"))

    # Randomise which 20 nodes are initial infectives — independent of the
    # network draw, so seed-position bias doesn't line up with hub position.
    set.seed(init_seed)
    seeded <- sample(seq_len(N), init_I)
    set.seed(NULL)
    I_ini_vec <- integer(N); I_ini_vec[seeded] <- 1L

    message(glue("Running SIR: pa = {pa}, network_seed = {nseed} ..."))
    raw <- run_stoch_adj_sparse(
      adj,
      beta           = N * beta_net / mean_k,
      t              = t_horizon,
      I_ini          = I_ini_vec,
      groups         = list(all = seq_len(N)),
      susceptibility = rep(1, N),
      gamma          = gamma,
      dt             = dt,
      timepoints     = timepoints,
      n_sim          = n_rep,
      cores          = 8)
    setDT(raw)

    # Cumulative NEW infections (excludes the seeded infectives; the odin
    # model initialises C = 0 and increments only from n_SI).
    dat <- raw[, .(time, sim, cases = C_all)]
    dat[, pl_alpha     := pa]
    dat[, network_seed := nseed]
    dat[, cv2_deg      := cv2_deg]
    all_rows[[length(all_rows) + 1L]] <- dat
  }
}

results <- rbindlist(all_rows)
fwrite(results, file.path(out_dir, "results.csv"))

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

# Per-replicate time to reach `target_cases`: first t where cumulative
# cases cross the threshold. NA if the replicate never gets there.
time_to_target <- results[, {
  hit <- which(cases >= target_cases)
  list(t_hit = if (length(hit)) time[hit[1L]] else NA_real_)
}, by = .(pl_alpha, network_seed, sim)]
fwrite(time_to_target, file.path(out_dir, "time_to_target.csv"))

# Median + 25-75% quantile band across replicates within each
# (pl_alpha, network_seed, time).
summary_dt <- results[, .(cases_med = median(cases),
                          cases_lo  = quantile(cases, 0.25),
                          cases_hi  = quantile(cases, 0.75)),
                      by = .(time, pl_alpha, network_seed)]

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

summary_dt[, pl_lbl := factor(sprintf("Pareto exp. = %s", format(pl_alpha, trim = TRUE)),
                              levels = sprintf("Pareto exp. = %s",
                                               format(pl_alphas, trim = TRUE)))]
time_to_target[, pl_lbl := factor(sprintf("Pareto exp. = %s", format(pl_alpha, trim = TRUE)),
                                   levels = sprintf("Pareto exp. = %s",
                                                    format(pl_alphas, trim = TRUE)))]

pal <- brewer.pal(max(3L, n_networks), "Dark2")[seq_len(n_networks)]

# 1. Cumulative cases vs time, one facet per pl_alpha, one coloured line
#    per network realisation with 25-75% ribbon across dust replicates.
p_cum <- ggplot(summary_dt,
                aes(x = time, y = cases_med,
                    colour = factor(network_seed),
                    fill   = factor(network_seed),
                    group  = network_seed)) +
  geom_ribbon(aes(ymin = cases_lo, ymax = cases_hi),
              alpha = 0.18, colour = NA) +
  geom_line(size = 1) +
  geom_hline(yintercept = target_cases, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ pl_lbl, scales = "free_y") +
  scale_colour_manual(name = "network realisation", values = pal) +
  scale_fill_manual(name   = "network realisation", values = pal) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(x = "t",
       y = "Cumulative new cases",
       title = glue("Early epidemic — N = {N}, R0 = {R0}, init_I = {init_I}, ",
                    "n_rep = {n_rep} per (pa, network)"))
ggsave(file.path(out_dir, "cumulative_cases.png"), p_cum,
       width = 11, height = 8, dpi = 130)

# 2. Same as above but log y-axis — exponential growth shows as straight
#    lines, so slopes across networks are visually comparable.
p_cum_log <- p_cum +
  scale_y_log10(labels = scales::comma) +
  labs(y = "Cumulative new cases (log scale)")
ggsave(file.path(out_dir, "cumulative_cases_log.png"), p_cum_log,
       width = 11, height = 8, dpi = 130)

# 3. Time to reach `target_cases`, dodged by network realisation within
#    each pl_alpha. Jitter overlays individual replicates.
p_hit <- ggplot(time_to_target[!is.na(t_hit)],
                aes(x = factor(pl_alpha), y = t_hit,
                    colour = factor(network_seed))) +
  geom_boxplot(outlier.shape = NA,
               position = position_dodge(width = 0.7),
               width = 0.55) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.12,
                                              dodge.width = 0.7),
              alpha = 0.3, size = 0.6) +
  scale_colour_manual(name = "network realisation", values = pal) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  labs(x = "Pareto exponent",
       y = glue("Time to reach {target_cases} cumulative cases"),
       title = glue("Time to threshold across networks ({n_networks} per pa, ",
                    "{n_rep} reps each)"))
ggsave(file.path(out_dir, "time_to_target.png"), p_hit,
       width = 10, height = 6, dpi = 130)

message(glue("Done. Outputs in {out_dir}/"))
