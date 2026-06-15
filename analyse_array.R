# Read all slurm-array RDS files and produce:
#   - debug table of effective sample sizes per MCMC run
#   - forest plots of VE and alpha:
#       * "separate":  each network contact matrix as its own row
#                      (10 allocations pooled within each network)
#       * "combined":  all networks pooled into one row
#     Frailty always pools all 10 allocations into one row.
#     Linear / sir are single rows each.
#
# Usage: Rscript analyse_array.R [results_dir] [out_dir]

library(data.table)
library(ggplot2)
library(dplyr)
library(coda)
library(glue)

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "output/array_results"
out_dir     <- if (length(args) >= 2) args[2] else file.path(results_dir, "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Load and flatten
# ---------------------------------------------------------------------------

files <- list.files(results_dir, pattern = "^results_.*\\.RDS$", full.names = TRUE)
if (length(files) == 0) stop("No results_*.RDS files in ", results_dir)
message(glue("Loading {length(files)} RDS files from {results_dir}"))

flat <- unlist(lapply(files, readRDS), recursive = FALSE)
ok   <- Filter(function(x) !is.null(x$result), flat)
message(glue("{length(ok)}/{length(flat)} configs returned a result"))

# ---------------------------------------------------------------------------
# Per-config extraction: VE samples at final time, alpha posterior, ESS
# ---------------------------------------------------------------------------

extract <- function(r, method_pref = "full") {
  cfg  <- r$config
  pars <- r$result$params
  ve   <- r$result$VE
  setDT(ve)

  t_final <- cfg$t

  # Normalise to a single "VE" column at the final time. Two output shapes:
  # - linear wide:   columns t, C1, C2, full, CRR, sim
  # - long format:   columns t, eate, method, sim (and maybe num/denom)
  if ("method" %in% names(ve)) {
    sel <- ve[t == t_final & method == method_pref, ]
    if (!nrow(sel)) sel <- ve[t == t_final & method == "CRR", ]
    ve_samples <- 1 - sel$eate
  } else {
    sel <- ve[t == t_final, ]
    col <- if (method_pref %in% names(sel)) method_pref else "CRR"
    ve_samples <- 1 - sel[[col]]
  }

  alpha_samples <- as.numeric(pars[, 2])

  list(
    name            = r$name,
    model_type      = cfg$model_type %||% "legacy",
    network_seed    = cfg$network_seed    %||% NA_integer_,
    allocation_seed = cfg$allocation_seed %||% NA_integer_,
    ve_samples      = ve_samples,
    alpha_samples   = alpha_samples,
    n_post          = nrow(pars),
    n_ve            = length(ve_samples),
    ess_beta        = as.numeric(coda::effectiveSize(pars[, 1])),
    ess_alpha       = as.numeric(coda::effectiveSize(pars[, 2]))
  )
}

infos <- lapply(ok, extract)

# ---------------------------------------------------------------------------
# Debug: ESS table
# ---------------------------------------------------------------------------

ess_table <- rbindlist(lapply(infos, function(x) {
  data.table(name            = x$name,
             model_type      = x$model_type,
             network_seed    = x$network_seed,
             allocation_seed = x$allocation_seed,
             n_post          = x$n_post,
             n_ve            = x$n_ve,
             ess_beta        = round(x$ess_beta, 1),
             ess_alpha       = round(x$ess_alpha, 1))
}))
setorder(ess_table, model_type, network_seed, allocation_seed)

ess_path <- file.path(out_dir, "ess_table.csv")
fwrite(ess_table, ess_path)
message(glue("Wrote {ess_path}"))

cat("\n=== Effective sample sizes (per MCMC run) ===\n")
print(ess_table)
cat("\nSummary of ESS across runs:\n")
print(ess_table[, .(min_beta  = min(ess_beta),
                    med_beta  = median(ess_beta),
                    max_beta  = max(ess_beta),
                    min_alpha = min(ess_alpha),
                    med_alpha = median(ess_alpha),
                    max_alpha = max(ess_alpha)),
                by = model_type])

# ---------------------------------------------------------------------------
# Pool samples by group and summarise (median, 95% interval)
# ---------------------------------------------------------------------------

# Group label functions: "separate" gives one row per network contact matrix
# (allocations pooled within); "combined" pools all networks into one row.
label_separate <- function(x) {
  switch(x$model_type,
         network = sprintf("network_n%02d", x$network_seed),
         frailty = "frailty",
         legacy  = x$name)
}
label_combined <- function(x) {
  switch(x$model_type,
         network = "network_all",
         frailty = "frailty",
         legacy  = x$name)
}

# Ordering for the y-axis: legacy first, then frailty, then networks
order_key <- function(label) {
  if (label %in% c("linear", "sir"))    return(paste0("0_", label))
  if (label == "frailty")               return("1_frailty")
  if (label == "network_all")           return("2_network_all")
  paste0("3_", label)
}

summarise_by_group <- function(infos, group_fn, field, q = c(0.025, 0.5, 0.975)) {
  pooled <- list()
  for (x in infos) {
    g <- group_fn(x)
    pooled[[g]] <- c(pooled[[g]], x[[field]])
  }
  rbindlist(lapply(names(pooled), function(g) {
    s <- pooled[[g]]
    s <- s[is.finite(s)]
    if (!length(s)) return(NULL)
    qs <- quantile(s, q)
    data.table(group  = g,
               n      = length(s),
               q025   = qs[[1]],
               median = qs[[2]],
               q975   = qs[[3]])
  }))[order(sapply(group, order_key))]
}

ve_sep    <- summarise_by_group(infos, label_separate, "ve_samples")
ve_comb   <- summarise_by_group(infos, label_combined, "ve_samples")
alpha_sep <- summarise_by_group(infos, label_separate, "alpha_samples")
alpha_comb<- summarise_by_group(infos, label_combined, "alpha_samples")

# ---------------------------------------------------------------------------
# Forest plots
# ---------------------------------------------------------------------------

forest_plot <- function(df, title, xlab, vline = NULL) {
  df[, group := factor(group, levels = rev(df$group))]
  p <- ggplot(df, aes(y = group, x = median)) +
    geom_point(size = 2.6) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.2) +
    geom_text(aes(label = sprintf("n=%d", n)),
              x = Inf, hjust = 1.1, size = 3, colour = "grey40") +
    theme_minimal(base_size = 13) +
    labs(x = xlab, y = NULL, title = title)
  if (!is.null(vline)) p <- p + geom_vline(xintercept = vline, linetype = "dashed", colour = "grey50")
  p
}

ggsave(file.path(out_dir, "forest_VE_separate.png"),
       forest_plot(ve_sep,  "VE — networks shown per contact matrix (allocations pooled)",
                   "VE = 1 - EATE", vline = 0),
       width = 8, height = 7, dpi = 130)

ggsave(file.path(out_dir, "forest_VE_combined.png"),
       forest_plot(ve_comb, "VE — all networks pooled",
                   "VE = 1 - EATE", vline = 0),
       width = 8, height = 3.5, dpi = 130)

ggsave(file.path(out_dir, "forest_alpha_separate.png"),
       forest_plot(alpha_sep,  "alpha (susceptibility[2]) — networks per contact matrix",
                   "alpha"),
       width = 8, height = 7, dpi = 130)

ggsave(file.path(out_dir, "forest_alpha_combined.png"),
       forest_plot(alpha_comb, "alpha (susceptibility[2]) — networks pooled",
                   "alpha"),
       width = 8, height = 3.5, dpi = 130)

message(glue("Wrote forest plots to {out_dir}/"))
message("Done.")
