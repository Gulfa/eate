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
    name            = as.character(r$name),
    model_type      = as.character(r$model_type),
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
# Group labels — three consistent aggregation levels:
#   L1 all_splits   : one row per job (every allocation, every network)
#   L2 pool_allocs  : one row per (model | network_seed), allocations pooled
#   L3 pool_nets    : one row per model, all allocations AND networks pooled
# ---------------------------------------------------------------------------

.frailty_short <- function(mt) sub("^sir_", "", mt)   # sir_sus_frailty -> sus_frailty

labels_L1_all_splits <- function(r) {
  switch(r$model_type,
         network           = sprintf("network_n%02d_a%02d", r$network_seed, r$allocation_seed),
         sir_sus_frailty   = sprintf("sus_frailty_a%02d",   r$allocation_seed),
         sir_trans_frailty = sprintf("trans_frailty_a%02d", r$allocation_seed),
         r$model_type)
}
labels_L2_pool_allocs <- function(r) {
  switch(r$model_type,
         network           = sprintf("network_n%02d", r$network_seed),
         sir_sus_frailty   = "sus_frailty",
         sir_trans_frailty = "trans_frailty",
         r$model_type)
}
labels_L3_pool_nets <- function(r) {
  switch(r$model_type,
         network           = "network_all",
         sir_sus_frailty   = "sus_frailty",
         sir_trans_frailty = "trans_frailty",
         r$model_type)
}

order_key <- function(label) {
  if (label %in% c("linear", "sir"))    return(paste0("0_", label))
  if (grepl("frailty", label))          return(paste0("1_", label))
  if (label == "network_all")           return("2_network_all")
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

# ---------------------------------------------------------------------------
# VE at t* summariser (uses ve_uncertainty draws for CIs)
# ---------------------------------------------------------------------------

# Determine t*
t_star_ve <- max(unlist(lapply(ok, function(r) {
  if (!is.null(r$ve_uncertainty) && nrow(r$ve_uncertainty)) max(r$ve_uncertainty$t) else NA_real_
})), na.rm = TRUE)
if (!is.finite(t_star_ve)) t_star_ve <- max(unlist(lapply(ok, function(r) max(r$ve$t, na.rm = TRUE))), na.rm = TRUE)
message(glue("Using t* = {t_star_ve} for VE forest plots"))

summarise_ve_by <- function(ok, group_fn, t_target) {
  draws <- rbindlist(lapply(ok, function(r) {
    if (is.null(r$ve_uncertainty) || !nrow(r$ve_uncertainty)) return(NULL)
    v <- r$ve_uncertainty[method == "full_stoch" & t == t_target]
    if (!nrow(v)) return(NULL)
    data.table(group = group_fn(r), VE = 1 - v$eate)
  }))
  if (!nrow(draws)) return(data.table())
  s <- draws[, .(n = .N,
                 estimate = median(VE, na.rm = TRUE),
                 lo       = quantile(VE, 0.025, na.rm = TRUE),
                 hi       = quantile(VE, 0.975, na.rm = TRUE)),
             by = group]
  s[order(sapply(group, order_key))]
}

# ---------------------------------------------------------------------------
# Level-driven forest plots for beta, alpha, VE at t*
# ---------------------------------------------------------------------------

levels_def <- list(
  all_splits  = labels_L1_all_splits,
  pool_allocs = labels_L2_pool_allocs,
  pool_nets   = labels_L3_pool_nets
)
level_titles <- list(
  all_splits  = "all splits (per job)",
  pool_allocs = "allocations pooled (per network / model)",
  pool_nets   = "everything pooled (per model)"
)

for (lvl in names(levels_def)) {
  lbl_fn <- levels_def[[lvl]]

  df_beta  <- summarise_param(ok, lbl_fn, "beta")
  df_alpha <- summarise_param(ok, lbl_fn, "alpha")
  df_ve    <- summarise_ve_by(ok, lbl_fn, t_star_ve)

  # Scale height with row count so labels stay legible
  h <- max(3.5, 0.28 * nrow(df_beta))
  small <- if (lvl == "all_splits") theme(axis.text.y = element_text(size = 7)) else NULL

  ggsave(file.path(out_dir, glue("forest_beta_{lvl}.png")),
         forest_plot(df_beta, glue("beta — {level_titles[[lvl]]}"), "beta") + small,
         width = 8, height = h, dpi = 130, limitsize = FALSE)
  ggsave(file.path(out_dir, glue("forest_alpha_{lvl}.png")),
         forest_plot(df_alpha, glue("alpha — {level_titles[[lvl]]}"), "alpha", vline = 1) + small,
         width = 8, height = h, dpi = 130, limitsize = FALSE)
  if (nrow(df_ve)) {
    ggsave(file.path(out_dir, glue("forest_VE_t{t_star_ve}_{lvl}.png")),
           forest_plot(df_ve, glue("VE at t = {t_star_ve} — {level_titles[[lvl]]}"),
                       "VE = 1 - EATE", vline = 0) + small,
           width = 8, height = h, dpi = 130, limitsize = FALSE)
  }

  # Also save the underlying summary CSVs
  fwrite(df_beta,  file.path(out_dir, glue("forest_beta_{lvl}.csv")))
  fwrite(df_alpha, file.path(out_dir, glue("forest_alpha_{lvl}.csv")))
  if (nrow(df_ve)) fwrite(df_ve, file.path(out_dir, glue("forest_VE_t{t_star_ve}_{lvl}.csv")))
}

# ---------------------------------------------------------------------------
# VE trajectories
# ---------------------------------------------------------------------------

# Long VE table: one row per (job, t, sim, method)
ve_long <- rbindlist(lapply(ok, function(r) {
  v <- copy(r$ve)
  v[, name            := as.character(r$name)]
  v[, model_type      := as.character(r$model_type)]
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

# ---------------------------------------------------------------------------
# VE with propagated parameter uncertainty
# ---------------------------------------------------------------------------

# Long table: one row per (job, param_sample, allocation/sim, t, method)
ve_unc_long <- rbindlist(lapply(ok, function(r) {
  if (is.null(r$ve_uncertainty) || !nrow(r$ve_uncertainty)) return(NULL)
  v <- copy(r$ve_uncertainty)
  v[, name            := as.character(r$name)]
  v[, model_type      := as.character(r$model_type)]
  v[, network_seed    := r$network_seed    %||% NA_integer_]
  v[, allocation_seed := r$allocation_seed %||% NA_integer_]
  v
}), fill = TRUE)

if (nrow(ve_unc_long) > 0) {
  ve_unc_long[, VE := 1 - eate]
  ve_unc <- ve_unc_long[method == "full_stoch"]

  # Variance decomposition: we want a clean "parameter-only" band that
  # doesn't accidentally mix in outer-allocation spread.
  #   1. For each (job, param_sample, t): mean VE over inner sims.
  #   2. For each (job, t): take the 2.5/50/97.5 quantile across
  #      param_samples -> per-job parameter-only spread.
  #   3. Pool across jobs by averaging the within-job bounds (so the
  #      band is "typical per-allocation parameter uncertainty").
  # The allocation-only band is symmetric (collapse over param_sample
  # first, then quantile over allocations within model).
  # Total band: pool every (param_sample, sim, outer alloc) per (t, model).
  vu_jps <- ve_unc[, .(VE = mean(VE, na.rm = TRUE)),
                  by = .(t, model_type, network_seed, allocation_seed, param_sample)]

  param_per_job <- vu_jps[, .(med = median(VE),
                              lo  = quantile(VE, 0.025, na.rm = TRUE),
                              hi  = quantile(VE, 0.975, na.rm = TRUE)),
                          by = .(t, model_type, network_seed, allocation_seed)]
  param_band <- param_per_job[, .(VE_med = median(med),
                                  lo = mean(lo), hi = mean(hi),
                                  n  = .N),
                              by = .(t, model_type)]

  vu_alloc <- ve_unc[, .(VE = mean(VE, na.rm = TRUE)),
                    by = .(t, model_type, network_seed, allocation_seed, sim)]
  # Per-param allocation-only spread (mirror of the param logic but
  # collapsing over param_sample first).
  vu_aps <- ve_unc[, .(VE = mean(VE, na.rm = TRUE)),
                  by = .(t, model_type, network_seed, allocation_seed, sim)]
  alloc_band <- vu_aps[, .(VE_med = median(VE),
                           lo = quantile(VE, 0.025, na.rm = TRUE),
                           hi = quantile(VE, 0.975, na.rm = TRUE),
                           n  = .N),
                       by = .(t, model_type)]

  total_band <- ve_unc[, .(VE_med = median(VE),
                           lo = quantile(VE, 0.025, na.rm = TRUE),
                           hi = quantile(VE, 0.975, na.rm = TRUE),
                           n  = .N),
                       by = .(t, model_type)]

  param_band[, source := "param"]
  alloc_band[, source := "allocation"]
  total_band[, source := "total"]
  bands <- rbindlist(list(param_band, alloc_band, total_band))

  fwrite(bands, file.path(out_dir, "ve_uncertainty_bands.csv"))

  # Double-ribbon plot: inner (param-only) + outer (total). Allocation
  # band shown by colour overlay if you want, but for the headline plot
  # param + total is most informative.
  p_unc <- ggplot(total_band, aes(x = t, y = VE_med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                fill = "grey60", alpha = 0.4) +
    geom_ribbon(data = param_band,
                aes(ymin = lo, ymax = hi),
                fill = "steelblue", alpha = 0.5) +
    geom_line(size = 1, colour = "black") +
    facet_wrap(~ model_type, scales = "free_y") +
    theme_minimal(base_size = 13) +
    labs(x = "t", y = "VE = 1 - EATE",
         title = "VE(t) with posterior uncertainty",
         subtitle = "blue band = parameter uncertainty only; grey band = total (params + allocations)")
  ggsave(file.path(out_dir, "ve_trajectory_with_uncertainty.png"),
         p_unc, width = 11, height = 7, dpi = 130)

  # Final-time VE per model with the two CIs
  ve_unc_final <- bands[t == max(t)]
  fwrite(dcast(ve_unc_final, model_type ~ source,
               value.var = c("VE_med", "lo", "hi", "n")),
         file.path(out_dir, "ve_final_intervals.csv"))

  message("\n=== VE at t* with uncertainty (param vs total) ===")
  print(ve_unc_final[, .(model_type, source, VE_med, lo, hi, n)])
}

message(glue("\nWrote plots and CSVs to {out_dir}/"))
