# Read all run_fit_array.R outputs and produce:
#   - Forest plots of fitted beta and alpha per model, with the Gaussian
#     posterior central interval (mean +/- z * sd from J^-1 Sigma J^-T);
#     interval level set by `ci_level` below (default 80%).
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
library(RColorBrewer)

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "output/fit_array_results"
out_dir     <- if (length(args) >= 2) args[2] else file.path(results_dir, "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# Central-interval level for all posterior/uncertainty bands. 0.80 gives
# narrower, less jittery bands than 0.95: absolute wobble scales with z,
# so 80% cuts it ~35% (z 1.28 vs 1.96), and the endpoints sit in the bulk
# where the mean +/- z*SD Gaussian approximation holds best.
ci_level  <- 0.80
z_ci      <- qnorm(0.5 + ci_level / 2)     # 1.2816 at 80%, 1.96 at 95%
q_lo      <- (1 - ci_level) / 2            # 0.10 at 80%
q_hi      <- 1 - q_lo                      # 0.90 at 80%
ci_pct    <- round(100 * ci_level)         # for titles/labels

# ---------------------------------------------------------------------------
# Load + flatten
# ---------------------------------------------------------------------------

files <- list.files(results_dir, pattern = "^results_.*\\.RDS$", full.names = TRUE)
if (length(files) == 0) stop("No results_*.RDS files in ", results_dir)
message(glue("Loading {length(files)} RDS files from {results_dir}"))

flat <- unlist(lapply(files, readRDS), recursive = FALSE)
ok_all <- Filter(function(x) !is.null(x) && !is.null(x$fit), flat)
message(glue("{length(ok_all)}/{length(flat)} jobs returned a result"))

# Split by experiment_id (backwards-compat: NA -> "default")
experiment_ids <- vapply(ok_all,
                         function(r) r$experiment_id %||% "default",
                         character(1))
experiments_present <- unique(experiment_ids)
message(glue("Experiments found: {paste(experiments_present, collapse=', ')}"))

# ---------------------------------------------------------------------------
# Top-level helpers used by both the per-experiment loop and the
# post-loop combined plots. Keep OUTSIDE analyse_one_experiment so they
# remain visible after the loop.
# ---------------------------------------------------------------------------

.pa_str <- function(pa) sprintf("pa%s", format(pa, trim = TRUE))

labels_L1_all_splits <- function(r) {
  switch(r$model_type,
         network           = sprintf("network_%s_n%02d_a%02d",
                                     .pa_str(r$pl_alpha), r$network_seed, r$allocation_seed),
         sir_sus_frailty   = sprintf("sus_frailty_a%02d",   r$allocation_seed),
         sir_trans_frailty = sprintf("trans_frailty_a%02d", r$allocation_seed),
         r$model_type)
}
labels_L2_pool_allocs <- function(r) {
  switch(r$model_type,
         network           = sprintf("network_%s_n%02d",
                                     .pa_str(r$pl_alpha), r$network_seed),
         sir_sus_frailty   = "sus_frailty",
         sir_trans_frailty = "trans_frailty",
         r$model_type)
}
labels_L3_pool_nets <- function(r) {
  switch(r$model_type,
         network           = sprintf("network_%s_all", .pa_str(r$pl_alpha)),
         sir_sus_frailty   = "sus_frailty",
         sir_trans_frailty = "trans_frailty",
         r$model_type)
}

order_key <- function(label) {
  if (label %in% c("linear", "sir", "sir_separate")) return(paste0("0_", label))
  if (grepl("frailty", label))       return(paste0("1_", label))
  if (grepl("^network_pa.*_all$", label)) return(paste0("2_", label))
  paste0("3_", label)
}

# Human-readable labels for the group codes used by the forest plots.
# Handles L1 (job-specific), L2 (per network_seed), L3 (per pl_alpha)
# variants. Anything unrecognised falls through unchanged.
display_name <- function(x) {
  vapply(as.character(x), function(g) {
    if (g == "linear")       return("Linear")
    if (g == "sir")          return("SIR (homogeneous)")
    if (g == "sir_separate") return("SIR (separate pops)")
    if (grepl("^sus_frailty(_a[0-9]+)?$", g)) {
      s <- sub("^sus_frailty",  "SIR + sus. frailty",  g)
      return(sub("_a([0-9]+)$", " (alloc \\1)", s))
    }
    if (grepl("^trans_frailty(_a[0-9]+)?$", g)) {
      s <- sub("^trans_frailty", "SIR + trans. frailty", g)
      return(sub("_a([0-9]+)$", " (alloc \\1)", s))
    }
    m <- regmatches(g, regexec("^network_pa([0-9.]+)(_n([0-9]+))?(_a([0-9]+))?(_all)?$", g))[[1]]
    if (length(m) >= 2 && nzchar(m[2])) {
      pa   <- m[2]
      seed <- m[4]; alloc <- m[6]; all_tag <- m[7]
      base <- sprintf("Network (Pareto exp. = %s)", pa)
      if (nzchar(all_tag) || (!nzchar(seed) && !nzchar(alloc))) return(base)
      if (nzchar(seed) && nzchar(alloc))
        return(sprintf("%s, seed %s, alloc %s", base, seed, alloc))
      if (nzchar(seed))  return(sprintf("%s, seed %s", base, seed))
      if (nzchar(alloc)) return(sprintf("%s, alloc %s", base, alloc))
    }
    g
  }, character(1))
}

# ---------------------------------------------------------------------------
# Everything below runs per-experiment: takes a filtered ok list and an
# experiment-specific out_dir. Called once per experiment_id further down.
# ---------------------------------------------------------------------------

analyse_one_experiment <- function(ok, out_dir) {

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_dt <- rbindlist(lapply(ok, function(r) {
  data.table(
    name            = as.character(r$name),
    model_type      = as.character(r$model_type),
    network_seed    = r$network_seed    %||% NA_integer_,
    allocation_seed = r$allocation_seed %||% NA_integer_,
    beta            = r$fit$beta,
    alpha           = r$fit$alpha,
    loss            = r$fit$loss,
    loss_floor      = r$loss_floor %||% NA_real_,
    loss_chisq      = r$loss_chisq %||% NA_real_,
    convergence     = r$fit$convergence,
    sd_beta         = r$posterior_cov$sd[["beta"]],
    cor_ba          = r$posterior_cov$cov[1, 2] /
                      (r$posterior_cov$sd[["beta"]] *
                       r$posterior_cov$sd[["alpha"]]),
    sd_alpha        = r$posterior_cov$sd[["alpha"]],
    cov_ba          = r$posterior_cov$cov[1, 2]
  )
}))
fit_dt[, beta_lo  := beta  - z_ci * sd_beta]
fit_dt[, beta_hi  := beta  + z_ci * sd_beta]
fit_dt[, alpha_lo := alpha - z_ci * sd_alpha]
fit_dt[, alpha_hi := alpha + z_ci * sd_alpha]
# Fit-quality flags: loss relative to the MC noise floor, and a
# chi-square-scaled misfit test (~2 = perfect fit, > 6 = genuine misfit).
fit_dt[, loss_ratio := loss / loss_floor]
# NA (not FALSE) when chisq is unavailable (e.g. legacy results), so old
# fits aren't miscounted as misfits.
fit_dt[, acceptable := fifelse(is.finite(loss_chisq), loss_chisq <= 6, NA)]

fwrite(fit_dt, file.path(out_dir, "fit_summary.csv"))
message("\n=== Fit diagnostics ===")
print(fit_dt[, .(med_loss     = round(median(loss), 3),
                 med_floor    = round(median(loss_floor), 3),
                 med_chisq    = round(median(loss_chisq), 2),
                 max_chisq    = round(max(loss_chisq), 2),
                 n_misfit     = sum(!acceptable, na.rm = TRUE),
                 conv_nonzero = sum(convergence != 0),
                 n = .N),
             by = model_type])
n_bad <- sum(!fit_dt$acceptable, na.rm = TRUE)
if (n_bad > 0) {
  message(glue("\n{n_bad} fit(s) with loss_chisq > 6 (genuine misfit):"))
  print(fit_dt[!(acceptable), .(name, loss, loss_floor,
                                loss_chisq = round(loss_chisq, 2))][
                                order(-loss_chisq)][seq_len(min(.N, 20))])
}

# Group labels / order_key are defined at top-level scope (above the
# function definition) so they're visible both here and in the
# combined cross-experiment section that runs after the loop.

# ---------------------------------------------------------------------------
# Posterior-aware pooling by BUCKETING draws: pool the posterior draws
# (beta_k / alpha_k) stored in ve_uncertainty and take their marginal SD
# directly, so the beta/alpha forests use the EXACT same draws as the
# by-allocation histograms (and the VE forest, which already buckets its
# draws). This is the empirical law-of-total-variance: sd_total = sd(all
# draws) = sqrt(var_within + var_between). The two components are still
# reported for diagnostics (they need not reconcile exactly with sd_total,
# since the bucket weights by draw count while the decomposition weights
# jobs equally). Falls back to the point estimate for legacy results that
# lack stored draws.
# ---------------------------------------------------------------------------

summarise_param <- function(ok, group_fn, param) {
  col <- paste0(param, "_k")
  draws <- rbindlist(lapply(seq_along(ok), function(i) {
    r   <- ok[[i]]
    grp <- group_fn(r)
    vu  <- r$ve_uncertainty
    if (is.null(vu) || !nrow(vu) || !(col %in% names(vu)))
      return(data.table(group = grp, job = i, value = r$fit[[param]]))
    v  <- vu[method == "full_stoch"]
    if (!nrow(v)) v <- vu
    dd <- unique(v[, .(param_sample, value = get(col))])
    data.table(group = grp, job = i, value = dd$value)
  }))

  # Within/between decomposition (equal weight per job), for diagnostics.
  per_job <- draws[, .(job_mean = mean(value, na.rm = TRUE),
                       job_var  = var(value,  na.rm = TRUE)),
                   by = .(group, job)]
  decomp <- per_job[, .(sd_within  = sqrt(mean(job_var, na.rm = TRUE)),
                        sd_between = if (.N > 1) sd(job_mean) else 0,
                        n_jobs     = .N),
                    by = group]
  # Bucketed marginal SD (matches the histogram exactly).
  s <- draws[, .(n        = .N,
                 estimate = mean(value, na.rm = TRUE),
                 sd_total = sd(value,   na.rm = TRUE)),
             by = group]
  s <- merge(s, decomp, by = "group")
  s[, lo := estimate - z_ci * sd_total]
  s[, hi := estimate + z_ci * sd_total]
  s[order(sapply(group, order_key))]
}

forest_plot <- function(df, title, xlab, vline = NULL) {
  ord   <- rev(as.character(df$group))
  labs  <- setNames(display_name(ord), ord)
  df[, group := factor(as.character(group), levels = ord, labels = labs[ord])]
  p <- ggplot(df, aes(y = group, x = estimate)) +
    geom_point(size = 2.6) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA)) +
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
  # SD-based CI: mean +/- z_ci * SD(pooled). Consistent with the Laplace/
  # MVN posterior and stabler at small K_post_samples than empirical
  # tail quantiles (noise scales ~ 1/sqrt(2K) instead of ~ 1/sqrt(K)).
  s <- draws[, .(n        = .N,
                 estimate = mean(VE, na.rm = TRUE),
                 lo       = mean(VE, na.rm = TRUE) - z_ci * sd(VE, na.rm = TRUE),
                 hi       = mean(VE, na.rm = TRUE) + z_ci * sd(VE, na.rm = TRUE)),
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
         forest_plot(df_beta, glue("beta — {level_titles[[lvl]]} ({ci_pct}% CI)"), "beta") + small,
         width = 8, height = h, dpi = 130, limitsize = FALSE)
  ggsave(file.path(out_dir, glue("forest_alpha_{lvl}.png")),
         forest_plot(df_alpha, glue("alpha — {level_titles[[lvl]]} ({ci_pct}% CI)"), "alpha") + small,
         width = 8, height = h, dpi = 130, limitsize = FALSE)
  if (nrow(df_ve)) {
    ggsave(file.path(out_dir, glue("forest_VE_t{t_star_ve}_{lvl}.png")),
           forest_plot(df_ve, glue("VE at t = {t_star_ve} — {level_titles[[lvl]]} ({ci_pct}% CI)"),
                       "VE = 1 - EATE") + small,
           width = 8, height = h, dpi = 130, limitsize = FALSE)
  }

  # Also save the underlying summary CSVs
  fwrite(df_beta,  file.path(out_dir, glue("forest_beta_{lvl}.csv")))
  fwrite(df_alpha, file.path(out_dir, glue("forest_alpha_{lvl}.csv")))
  if (nrow(df_ve)) fwrite(df_ve, file.path(out_dir, glue("forest_VE_t{t_star_ve}_{lvl}.csv")))

  # Two-column combined VE | alpha at this level. Shared y-axis of
  # model groups, free x-scales so VE and alpha get their own ranges.
  # No plot title.
  if (nrow(df_ve) && nrow(df_alpha)) {
    ve_side    <- copy(df_ve)[,    metric := "VE"]
    alpha_side <- copy(df_alpha)[, metric := "alpha"]
    keep <- c("group", "n", "estimate", "lo", "hi", "metric")
    side <- rbind(ve_side[, ..keep], alpha_side[, ..keep])
    lvl_groups <- unique(as.character(side$group))
    lvl_groups <- lvl_groups[order(sapply(lvl_groups, order_key))]
    disp_labs  <- setNames(display_name(lvl_groups), lvl_groups)
    side[, group := factor(as.character(group),
                            levels = rev(lvl_groups),
                            labels = disp_labs[rev(lvl_groups)])]
    side[, metric := factor(metric, levels = c("VE", "alpha"))]

    p_side <- ggplot(side, aes(y = group, x = estimate)) +
      geom_point(size = 2.6) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
      facet_wrap(~ metric, scales = "free_x") +
      theme_bw(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey95", colour = NA)) +
      labs(x = NULL, y = NULL)
    ggsave(file.path(out_dir, glue("forest_VE_alpha_{lvl}.png")),
           p_side, width = 11, height = h, dpi = 130, limitsize = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Histograms of the K posterior draws of VE / alpha / beta, filled by
# allocation. For each job, ve_uncertainty holds K posterior samples
# (param_sample), each tagged with its drawn beta_k / alpha_k; VE is the
# mean of 1 - eate over that draw's inner sims. Stacking the K draws by
# allocation shows how each allocation contributes to the pooled VE /
# alpha / beta distribution. Panels are per (model, metric) so each
# parameter keeps its own scale.
# ---------------------------------------------------------------------------

draws_dt <- rbindlist(lapply(seq_along(ok), function(i) {
  r <- ok[[i]]
  if (is.null(r$ve_uncertainty) || !nrow(r$ve_uncertainty)) return(NULL)
  v <- r$ve_uncertainty[method == "full_stoch" & t == t_star_ve]
  if (!nrow(v)) return(NULL)
  has_bk <- all(c("beta_k", "alpha_k") %in% names(v))
  # One row per posterior draw: beta_k / alpha_k are constant within a
  # param_sample; VE averaged over the draw's inner sims -> K per job.
  agg <- v[, .(VE    = mean(1 - eate, na.rm = TRUE),
               beta  = if (has_bk) beta_k[1]  else r$fit$beta,
               alpha = if (has_bk) alpha_k[1] else r$fit$alpha),
           by = param_sample]
  data.table(job             = i,
             model_type      = as.character(r$model_type),
             allocation_seed = r$allocation_seed %||% NA_integer_,
             network_seed    = r$network_seed    %||% NA_integer_,
             pl_alpha        = r$pl_alpha        %||% NA_real_,
             VE = agg$VE, alpha = agg$alpha, beta = agg$beta)
}))

if (!nrow(draws_dt)) {
  message("No ve_uncertainty draws available; skipping allocation histograms.")
} else {
  hist_dt <- melt(draws_dt,
                  id.vars       = c("model_type", "allocation_seed",
                                    "network_seed", "pl_alpha"),
                  measure.vars  = c("VE", "alpha", "beta"),
                  variable.name = "metric", value.name = "value")
  hist_dt[, metric := factor(metric, levels = c("VE", "alpha", "beta"))]

  dark2_pal <- function(n) if (n <= 8L) brewer.pal(max(n, 3L), "Dark2")[seq_len(n)]
                           else colorRampPalette(brewer.pal(8L, "Dark2"))(n)

  # Stacked histogram of the K posterior draws, faceted by `rowvar` x metric
  # (free scales) and filled by `fill_col`. Shows how each configuration
  # contributes to the pooled VE / alpha / beta distribution.
  save_hist <- function(dt, rowvar, fill_col, fill_name, title, file) {
    if (!nrow(dt)) return(invisible(NULL))
    d <- copy(dt)
    d[, rowf  := factor(get(rowvar))]
    d[, fillf := factor(get(fill_col))]
    pal <- dark2_pal(length(levels(d$fillf)))

    # Central 80% (10/90) and 95% (2.5/97.5) percentile edges of the pooled
    # draws, per panel (rowf x metric), as vertical lines.
    qs <- d[, {
      p <- quantile(value, c(0.025, 0.10, 0.90, 0.975), na.rm = TRUE)
      list(interval = c("95%", "80%", "80%", "95%"), x = as.numeric(p))
    }, by = .(rowf, metric)]
    qs[, interval := factor(interval, levels = c("80%", "95%"))]

    p <- ggplot(d, aes(x = value, fill = fillf)) +
      geom_histogram(bins = 40, position = "stack", colour = "white", linewidth = 0.1) +
      geom_vline(data = qs, aes(xintercept = x, linetype = interval),
                 inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
      facet_wrap(vars(rowf, metric), scales = "free", ncol = 3) +
      scale_fill_manual(name = fill_name, values = pal, na.value = "grey60") +
      scale_linetype_manual(name = "central interval",
                            values = c("80%" = "dashed", "95%" = "dotted")) +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey95", colour = NA)) +
      labs(x = NULL, y = "count", title = title)
    ggsave(file, p, width = 11,
           height = max(4, 2.6 * length(levels(d$rowf))), dpi = 130,
           limitsize = FALSE)
  }

  # Non-network models: contribution of each allocation.
  save_hist(hist_dt[model_type != "network"],
            rowvar = "model_type", fill_col = "allocation_seed",
            fill_name = "allocation",
            title = "Posterior draws (K per allocation) of VE / alpha / beta",
            file = file.path(out_dir, "hist_VE_alpha_beta_by_allocation.png"))

  # Network: contribution of each network configuration (contact matrix),
  # faceted by Pareto exponent since pl_alpha is a separate structural axis.
  save_hist(hist_dt[model_type == "network"],
            rowvar = "pl_alpha", fill_col = "network_seed",
            fill_name = "network_seed",
            title = "Network: posterior draws of VE / alpha / beta by network config (row = Pareto exp.)",
            file = file.path(out_dir, "hist_VE_alpha_beta_network_by_config.png"))

  fwrite(hist_dt, file.path(out_dir, "hist_VE_alpha_beta_draws.csv"))

  # -------------------------------------------------------------------------
  # Between-allocation spread per model: sd of the per-allocation MEAN of
  # each metric, plus how tightly VE tracks alpha across allocations.
  # Explains "alpha varies more across allocations than VE": VE is anchored
  # by the fitted data, alpha is latent, so sd_between(VE) << sd_between(alpha)
  # with a shallow dVE/dalpha slope. Grouped by (model, pl_alpha) so network
  # splits by Pareto exponent; for non-network models pl_alpha is NA and the
  # spread is genuinely across allocations.
  # -------------------------------------------------------------------------
  job_means <- draws_dt[, .(VE    = mean(VE,    na.rm = TRUE),
                            alpha = mean(alpha, na.rm = TRUE),
                            beta  = mean(beta,  na.rm = TRUE)),
                        by = .(model_type, pl_alpha, job)]
  between_tbl <- job_means[, {
    has_var <- .N > 1
    fit_ok  <- .N > 2 && sd(alpha) > 0 && sd(VE) > 0
    list(n_alloc          = .N,
         sd_between_beta  = if (has_var) sd(beta)  else 0,
         sd_between_alpha = if (has_var) sd(alpha) else 0,
         sd_between_VE    = if (has_var) sd(VE)    else 0,
         cor_alpha_VE     = if (fit_ok) cor(alpha, VE) else NA_real_,
         dVE_dalpha       = if (fit_ok) unname(coef(lm(VE ~ alpha))[2]) else NA_real_)
  }, by = .(model_type, pl_alpha)][order(model_type, pl_alpha)]
  fwrite(between_tbl, file.path(out_dir, "between_allocation_spread.csv"))
  message("\n=== Between-allocation spread (sd of per-allocation means) ===")
  print(between_tbl[, lapply(.SD, function(x)
                             if (is.numeric(x)) round(x, 4) else x)])
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

# Per-model VE summary at the final time, with ci_pct% interval across allocations
ve_final <- ve_alloc[t == max(t),
                     .(VE_med   = median(VE),
                       VE_lo    = quantile(VE, q_lo),
                       VE_hi    = quantile(VE, q_hi),
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
  #   2. For each (job, t): mean +/- z_ci * SD across
  #      param_samples -> per-job parameter-only spread.
  #   3. Pool across jobs by averaging the within-job bounds (so the
  #      band is "typical per-allocation parameter uncertainty").
  # The allocation-only band is symmetric (collapse over param_sample
  # first, then quantile over allocations within model).
  # Total band: pool every (param_sample, sim, outer alloc) per (t, model).
  vu_jps <- ve_unc[, .(VE = mean(VE, na.rm = TRUE)),
                  by = .(t, model_type, network_seed, allocation_seed, param_sample)]

  # Parameter-only band: mean +/- z_ci * SD across param_samples (per job,
  # per t), then average bounds across jobs. SD-based CI is stabler than
  # empirical tail quantiles at small K_post_samples.
  param_per_job <- vu_jps[, .(med = mean(VE, na.rm = TRUE),
                              lo  = mean(VE, na.rm = TRUE) -
                                    z_ci * sd(VE, na.rm = TRUE),
                              hi  = mean(VE, na.rm = TRUE) +
                                    z_ci * sd(VE, na.rm = TRUE)),
                          by = .(t, model_type, network_seed, allocation_seed)]
  param_band <- param_per_job[, .(VE_med = mean(med),
                                  lo = mean(lo), hi = mean(hi),
                                  n  = .N),
                              by = .(t, model_type)]

  vu_alloc <- ve_unc[, .(VE = mean(VE, na.rm = TRUE)),
                    by = .(t, model_type, network_seed, allocation_seed, sim)]
  # Per-param allocation-only spread (mirror of the param logic but
  # collapsing over param_sample first). Kept as quantile-based because
  # this is spread across ALLOCATIONS — a structural axis, not a Laplace
  # posterior — and typically not Gaussian.
  vu_aps <- ve_unc[, .(VE = mean(VE, na.rm = TRUE)),
                  by = .(t, model_type, network_seed, allocation_seed, sim)]
  alloc_band <- vu_aps[, .(VE_med = median(VE),
                           lo = quantile(VE, q_lo, na.rm = TRUE),
                           hi = quantile(VE, q_hi, na.rm = TRUE),
                           n  = .N),
                       by = .(t, model_type)]

  # Total band: SD-based on the pooled param x allocation x sim draws.
  total_band <- ve_unc[, .(VE_med = mean(VE, na.rm = TRUE),
                           lo = mean(VE, na.rm = TRUE) - z_ci * sd(VE, na.rm = TRUE),
                           hi = mean(VE, na.rm = TRUE) + z_ci * sd(VE, na.rm = TRUE),
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
         title = glue("VE(t) with posterior uncertainty ({ci_pct}% intervals)"),
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

}   # end analyse_one_experiment

# ---------------------------------------------------------------------------
# Run per-experiment
# ---------------------------------------------------------------------------

for (exp_id in experiments_present) {
  ok_exp <- ok_all[experiment_ids == exp_id]
  message(glue("\n========== Experiment: {exp_id}  ({length(ok_exp)} jobs) =========="))
  exp_out_dir <- file.path(out_dir, exp_id)
  tryCatch(
    analyse_one_experiment(ok_exp, exp_out_dir),
    error = function(e) message("Error in experiment ", exp_id, ": ", conditionMessage(e))
  )
}
message("\nAll experiments done.")

# ---------------------------------------------------------------------------
# Combined cross-experiment L3 (pool_nets) plots
# ---------------------------------------------------------------------------
# One dodged forest per parameter (beta, alpha, VE at t*) with all
# experiments overlaid, coloured by experiment_id. Uses the per-
# experiment L3 CSVs already written inside analyse_one_experiment.

combined_dir <- file.path(out_dir, "_combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

# Read one experiment's L3 CSV for a parameter, tag with experiment_id.
.read_l3 <- function(eid, param) {
  csv <- file.path(out_dir, eid, glue("forest_{param}_pool_nets.csv"))
  if (!file.exists(csv)) return(NULL)
  dt <- fread(csv)
  dt[, experiment_id := eid]
  dt[]
}

# VE has t_star embedded in the filename (glob for it).
.read_l3_ve <- function(eid) {
  files <- list.files(file.path(out_dir, eid),
                      pattern = "^forest_VE_t.*_pool_nets\\.csv$",
                      full.names = TRUE)
  if (!length(files)) return(NULL)
  dt <- fread(files[1])
  t_val <- as.numeric(sub("^forest_VE_t([0-9.]+)_pool_nets\\.csv$", "\\1",
                          basename(files[1])))
  dt[, experiment_id := eid]
  dt[, t_star := t_val]
  dt[]
}

beta_combined  <- rbindlist(lapply(experiments_present, .read_l3, param = "beta"), fill = TRUE)
alpha_combined <- rbindlist(lapply(experiments_present, .read_l3, param = "alpha"), fill = TRUE)
ve_combined    <- rbindlist(lapply(experiments_present, .read_l3_ve), fill = TRUE)

# Consistent group ordering across experiments (use order_key from
# analyse_one_experiment scope — redefine here for the top-level plots).
.order_group <- function(dt) {
  if (!nrow(dt)) return(dt)
  lvl  <- unique(as.character(dt$group))
  lvl  <- lvl[order(sapply(lvl, order_key))]
  labs <- setNames(display_name(lvl), lvl)
  dt[, group := factor(as.character(group),
                        levels = rev(lvl), labels = labs[rev(lvl)])]
  dt
}
beta_combined  <- .order_group(beta_combined)
alpha_combined <- .order_group(alpha_combined)
ve_combined    <- .order_group(ve_combined)

# Dark2 palette by experiment.
combined_forest <- function(df, xlab, title, vline = NULL, out_file) {
  if (!nrow(df)) return(invisible(NULL))
  n_exp <- length(unique(df$experiment_id))
  pal   <- if (n_exp <= 8L) brewer.pal(max(n_exp, 3L), "Dark2")[seq_len(n_exp)]
           else colorRampPalette(brewer.pal(8L, "Dark2"))(n_exp)
  p <- ggplot(df, aes(y = group, x = estimate,
                       colour = experiment_id, group = experiment_id)) +
    geom_point(position = position_dodge(width = 0.6), size = 2.6) +
    geom_errorbarh(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = 0.6), height = 0.25) +
    scale_colour_manual(name = "experiment", values = pal) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA)) +
    labs(x = xlab, y = NULL, title = title)
  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed", colour = "grey50")
  h <- max(4, 0.55 * length(unique(df$group)))
  ggsave(out_file, p, width = 10, height = h, dpi = 130, limitsize = FALSE)
}

combined_forest(beta_combined,
                xlab = "beta",
                title = glue("beta (pool_nets) across experiments ({ci_pct}% CI)"),
                out_file = file.path(combined_dir, "combined_forest_beta.png"))
combined_forest(alpha_combined,
                xlab = "alpha",
                title = glue("alpha (pool_nets) across experiments ({ci_pct}% CI)"),
                out_file = file.path(combined_dir, "combined_forest_alpha.png"))
combined_forest(ve_combined,
                xlab = "VE = 1 - EATE",
                title = glue("VE at t* (pool_nets) across experiments ({ci_pct}% CI)"),
                out_file = file.path(combined_dir, "combined_forest_VE.png"))

fwrite(beta_combined,  file.path(combined_dir, "combined_beta.csv"))
fwrite(alpha_combined, file.path(combined_dir, "combined_alpha.csv"))
fwrite(ve_combined,    file.path(combined_dir, "combined_VE.csv"))

# Two-column combined plot: VE (left) + alpha (right), both
# cross-experiment, shared y-axis (model groups), colour by experiment.
if (nrow(ve_combined) && nrow(alpha_combined)) {
  ve_side    <- copy(ve_combined)[,    metric := "VE"]
  alpha_side <- copy(alpha_combined)[, metric := "alpha"]
  keep <- c("group", "n", "estimate", "lo", "hi", "experiment_id", "metric")
  side_by_side <- rbind(ve_side[, ..keep], alpha_side[, ..keep])

  # Shared y-axis: ve_combined and alpha_combined already carry
  # display-name factor levels in order_key order (via .order_group),
  # so we take the union and preserve first-seen order.
  lvl_ve    <- levels(ve_combined$group)
  lvl_alpha <- levels(alpha_combined$group)
  lvl       <- union(lvl_ve, lvl_alpha)
  side_by_side[, group := factor(as.character(group), levels = lvl)]
  side_by_side[, metric := factor(metric, levels = c("VE", "alpha"))]

  n_exp <- length(unique(side_by_side$experiment_id))
  pal   <- if (n_exp <= 8L) brewer.pal(max(n_exp, 3L), "Dark2")[seq_len(n_exp)]
           else colorRampPalette(brewer.pal(8L, "Dark2"))(n_exp)

  p_two <- ggplot(side_by_side, aes(y = group, x = estimate,
                                     colour = experiment_id,
                                     group = experiment_id)) +
    geom_point(position = position_dodge(width = 0.6), size = 2.6) +
    geom_errorbarh(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = 0.6), height = 0.25) +
    facet_wrap(~ metric, scales = "free_x") +
    scale_colour_manual(name = "experiment", values = pal) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA)) +
    labs(x = NULL, y = NULL,
         title = glue("VE (t*) and alpha across experiments (pool_nets, {ci_pct}% CI)"))

  h <- max(4, 0.55 * length(lvl))
  ggsave(file.path(combined_dir, "combined_forest_VE_alpha.png"),
         p_two, width = 12, height = h, dpi = 130, limitsize = FALSE)
}

message(glue("Wrote combined cross-experiment plots to {combined_dir}/"))
