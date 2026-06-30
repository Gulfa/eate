# Read all run_fit_array.R outputs and produce:
#   - Forest plots of fitted beta and alpha per model, with the Gaussian
#     posterior 95% CI (mean +/- 1.96 * sd from J^-1 Sigma J^-T).
#   - VE(t) trajectory plots, per model, method = "full_stoch", with
#     median + min/max ribbon across allocations / replicates.
#   - Network configs in two flavours: "separate" (per network_seed,
#     allocations pooled) and "combined" (all networks pooled).
#   - Fit-diagnostics CSV (loss, convergence, posterior sd).
#
# Usage: Rscript analyse_fit_array.R [results_dir] [out_dir]

library(data.table)
library(ggplot2)
library(dplyr)
library(glue)

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "output/fit_array_results"
out_dir     <- if (length(args) >= 2) args[2] else file.path(results_dir, "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Load + flatten
# ---------------------------------------------------------------------------

files <- list.files(results_dir, pattern = "^results_.*\\.RDS$", full.names = TRUE)
if (length(files) == 0) stop("No results_*.RDS files in ", results_dir)
message(glue("Loading {length(files)} RDS files from {results_dir}"))

flat <- unlist(lapply(files, readRDS), recursive = FALSE)
ok   <- Filter(function(x) !is.null(x) && !is.null(x$fit), flat)
message(glue("{length(ok)}/{length(flat)} jobs returned a result"))

# ---------------------------------------------------------------------------
# Per-config fit table
# ---------------------------------------------------------------------------

fit_dt <- rbindlist(lapply(ok, function(r) {
  data.table(
    name            = r$name,
    model_type      = r$model_type,
    network_seed    = r$network_seed    %||% NA_integer_,
    allocation_seed = r$allocation_seed %||% NA_integer_,
    beta            = r$fit$beta,
    alpha           = r$fit$alpha,
    loss            = r$fit$loss,
    convergence     = r$fit$convergence,
    sd_beta         = r$posterior_cov$sd[["beta"]],
    sd_alpha        = r$posterior_cov$sd[["alpha"]],
    cov_ba          = r$posterior_cov$cov[1, 2]
  )
}))
fit_dt[, beta_lo  := beta  - 1.96 * sd_beta]
fit_dt[, beta_hi  := beta  + 1.96 * sd_beta]
fit_dt[, alpha_lo := alpha - 1.96 * sd_alpha]
fit_dt[, alpha_hi := alpha + 1.96 * sd_alpha]

fwrite(fit_dt, file.path(out_dir, "fit_summary.csv"))
message("\n=== Fit diagnostics ===")
print(fit_dt[, .(min_loss = min(loss), med_loss = median(loss),
                 max_loss = max(loss),
                 conv_nonzero = sum(convergence != 0),
                 n = .N),
             by = model_type])

# ---------------------------------------------------------------------------
# Group labels for forest plots (mirrors analyse_array.R conventions).
# ---------------------------------------------------------------------------

label_separate <- function(r) {
  switch(r$model_type,
         network = sprintf("network_n%02d", r$network_seed),
         r$name)
}
label_combined <- function(r) {
  switch(r$model_type,
         network = "network_all",
         r$name)
}
label_per_alloc <- function(r) {
  switch(r$model_type,
         network = sprintf("network_n%02d_a%02d", r$network_seed, r$allocation_seed),
         r$name)
}

order_key <- function(label) {
  if (label %in% c("linear", "sir"))            return(paste0("0_", label))
  if (grepl("^sir_.*frailty$", label))          return(paste0("1_", label))
  if (label == "network_all")                   return("2_network_all")
  paste0("3_", label)
}

# ---------------------------------------------------------------------------
# Posterior-aware pooling: for a group with multiple fits (e.g. network
# allocations within a seed), report (mean of betas, sd combining across-
# allocation spread + per-fit posterior sd). When only one fit is in a
# group, just use that fit's posterior sd directly.
# ---------------------------------------------------------------------------

pool_param <- function(values, sds) {
  # Pooled mean
  mu <- mean(values)
  # Total variance = mean of within-fit variance + between-fit variance
  var_within  <- mean(sds^2, na.rm = TRUE)
  var_between <- if (length(values) > 1) var(values) else 0
  list(mean = mu, sd = sqrt(var_within + var_between), n = length(values))
}

summarise_param <- function(ok, group_fn, param) {
  groups <- list()
  for (r in ok) {
    g <- group_fn(r)
    groups[[g]] <- rbindlist(list(groups[[g]],
                                  data.table(value = r$fit[[param]],
                                             sd    = r$posterior_cov$sd[[param]])))
  }
  rbindlist(lapply(names(groups), function(g) {
    p <- pool_param(groups[[g]]$value, groups[[g]]$sd)
    data.table(group = g, n = p$n,
               estimate = p$mean,
               lo = p$mean - 1.96 * p$sd,
               hi = p$mean + 1.96 * p$sd)
  }))[order(sapply(group, order_key))]
}

forest_plot <- function(df, title, xlab, vline = NULL) {
  df[, group := factor(group, levels = rev(df$group))]
  p <- ggplot(df, aes(y = group, x = estimate)) +
    geom_point(size = 2.6) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_text(aes(label = sprintf("n=%d", n)),
              x = Inf, hjust = 1.1, size = 3, colour = "grey40") +
    theme_minimal(base_size = 13) +
    labs(x = xlab, y = NULL, title = title)
  if (!is.null(vline)) p <- p + geom_vline(xintercept = vline, linetype = "dashed", colour = "grey50")
  p
}

# --- beta ---
beta_sep  <- summarise_param(ok, label_separate, "beta")
beta_comb <- summarise_param(ok, label_combined, "beta")
ggsave(file.path(out_dir, "forest_beta_separate.png"),
       forest_plot(beta_sep,  "beta — per network (allocations pooled)", "beta"),
       width = 8, height = 7, dpi = 130)
ggsave(file.path(out_dir, "forest_beta_combined.png"),
       forest_plot(beta_comb, "beta — networks pooled", "beta"),
       width = 8, height = 3.5, dpi = 130)

# --- alpha ---
alpha_sep  <- summarise_param(ok, label_separate, "alpha")
alpha_comb <- summarise_param(ok, label_combined, "alpha")
ggsave(file.path(out_dir, "forest_alpha_separate.png"),
       forest_plot(alpha_sep,  "alpha (vac susceptibility) — per network",
                   "alpha", vline = 1),
       width = 8, height = 7, dpi = 130)
ggsave(file.path(out_dir, "forest_alpha_combined.png"),
       forest_plot(alpha_comb, "alpha (vac susceptibility) — networks pooled",
                   "alpha", vline = 1),
       width = 8, height = 3.5, dpi = 130)

# --- per-allocation alpha forest (tall plot, one row per MCMC fit) ---
alpha_per_alloc <- summarise_param(ok, label_per_alloc, "alpha")
n_rows <- nrow(alpha_per_alloc)
p_apa <- forest_plot(alpha_per_alloc,
                     "alpha — per allocation", "alpha", vline = 1) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(out_dir, "forest_alpha_per_allocation.png"),
       p_apa, width = 8, height = max(4, 0.15 * n_rows),
       dpi = 130, limitsize = FALSE)

# ---------------------------------------------------------------------------
# VE trajectories
# ---------------------------------------------------------------------------

# Long VE table: one row per (job, t, sim, method)
ve_long <- rbindlist(lapply(ok, function(r) {
  v <- copy(r$ve)
  v[, name            := r$name]
  v[, model_type      := r$model_type]
  v[, network_seed    := r$network_seed    %||% NA_integer_]
  v[, allocation_seed := r$allocation_seed %||% NA_integer_]
  v
}), fill = TRUE)
ve_long[, VE := 1 - eate]

# Aggregate to per-(model, network_seed, allocation_seed, t) means across sim
ve_alloc <- ve_long[method == "full_stoch",
                    .(VE = mean(VE, na.rm = TRUE)),
                    by = .(t, model_type, network_seed, allocation_seed)]

# "Separate" trajectory: median + min/max ribbon across allocations per
# (model, network_seed, t). Network gets one band per seed; other models
# have a single allocation so the ribbon collapses to a line.
ve_sep <- ve_alloc[, .(VE_med = median(VE),
                       VE_min = min(VE),
                       VE_max = max(VE),
                       n      = .N),
                   by = .(t, model_type, network_seed)]

# "Combined" trajectory: pool across all allocations (network across
# network_seed too) per (model, t).
ve_comb <- ve_alloc[, .(VE_med = median(VE),
                        VE_min = min(VE),
                        VE_max = max(VE),
                        n      = .N),
                    by = .(t, model_type)]

p_sep <- ggplot(ve_sep, aes(x = t, y = VE_med, group = network_seed,
                            colour = factor(network_seed))) +
  geom_ribbon(aes(ymin = VE_min, ymax = VE_max, fill = factor(network_seed)),
              alpha = 0.15, colour = NA) +
  geom_line(size = 1) +
  facet_wrap(~ model_type, scales = "free_y") +
  scale_colour_viridis_d(name = "network_seed", na.value = "black") +
  scale_fill_viridis_d(name   = "network_seed", na.value = "black") +
  theme_minimal(base_size = 13) +
  labs(x = "t", y = "VE = 1 - EATE (full_stoch)",
       title = "VE(t) per model — network shown per contact matrix (allocs pooled)")
ggsave(file.path(out_dir, "ve_trajectory_separate.png"),
       p_sep, width = 11, height = 7, dpi = 130)

p_comb <- ggplot(ve_comb, aes(x = t, y = VE_med, group = model_type,
                              colour = model_type, fill = model_type)) +
  geom_ribbon(aes(ymin = VE_min, ymax = VE_max), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal(base_size = 13) +
  labs(x = "t", y = "VE = 1 - EATE (full_stoch)",
       title = "VE(t) per model (allocations + networks pooled)")
ggsave(file.path(out_dir, "ve_trajectory_combined.png"),
       p_comb, width = 9, height = 5, dpi = 130)

# Per-model VE summary at the final time, with 95% interval across allocations
ve_final <- ve_alloc[t == max(t),
                     .(VE_med   = median(VE),
                       VE_q025  = quantile(VE, 0.025),
                       VE_q975  = quantile(VE, 0.975),
                       n        = .N),
                     by = .(model_type)]
fwrite(ve_final, file.path(out_dir, "ve_final.csv"))
message("\n=== VE at t* (allocations pooled) ===")
print(ve_final)

message(glue("\nWrote plots and CSVs to {out_dir}/"))
