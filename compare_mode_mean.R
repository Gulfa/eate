# compare_mode_mean.R
# ---------------------------------------------------------------------------
# For ONE parameter set (beta, alpha), run every model's simulator for n_sim
# stochastic realisations and compare the vaccinated/unvaccinated case ratio
# (C2/C1) at the MEAN vs at the MODE (the takeoff component).
#
# Why: the mean-matching fit sets E[C] = data; the kernel/mode fit sets
# mode(C) = data. They give the SAME parameters (and the same VE) only if the
# distribution's mode = its mean. This script shows, per model, whether that
# holds -- i.e. whether switching to kernel fitting would move the VE.
#
# Run (in the fhi container, from /fhi/eate):  Rscript compare_mode_mean.R
# ---------------------------------------------------------------------------

Sys.setenv(EATE_SOURCE_ONLY = "1")            # load defs without running the array
suppressWarnings(suppressMessages(source("run_fit_array.R")))
suppressMessages(library(data.table))

# ---- config ---------------------------------------------------------------
beta         <- 1.3
alpha        <- 0.5
n_sim        <- 5000
seed         <- 1
cores        <- 8
takeoff_frac <- 0.10        # total attack > this fraction of N counts as a "takeoff"

# Shared setup for all models (data_* unused here; we only simulate).
base <- modifyList(base_common, list(
  experiment_id = "diag",
  N_cont = 300, N_vac = 300,
  data_C1 = NA_real_, data_C2 = NA_real_,
  t_star = 10, I_ini_2g = c(1, 1), init_I_nw = 2,
  gamma = 1, dt = 0.05, inner_cores = cores))

# Per-model overrides (mirrors build_configs_for_experiment).
models <- list(
  linear            = list(model_type = "linear"),
  sir               = list(model_type = "sir"),
  sir_multisite     = list(model_type = "sir_multisite", n_sites = 4, site_icc = 0,
                           allocation_seed = 1),
  sir_sus_frailty   = list(model_type = "sir_sus_frailty", sd = 0.3, sd_trans = 0,
                           n_frailty = 10, frailty_amp = 2.5, allocation_seed = 1),
  sir_trans_frailty = list(model_type = "sir_trans_frailty", sd = 0, sd_trans = 0.3,
                           n_frailty = 10, frailty_amp = 2.5, allocation_seed = 1))
  # network (slow) commented out while testing sir_split_effect:
  # , network         = list(model_type = "network", pl_alpha = 1.5, mean_k = 6,
  #                          network_seed = 1, allocation_seed = 1))

# ---- run + summarise ------------------------------------------------------
N_tot <- base$N_cont + base$N_vac
raw   <- list()
summ  <- rbindlist(lapply(names(models), function(nm) {
  cfg <- materialise_cfg(modifyList(base, models[[nm]]))
  out <- build_simulator(cfg)(beta, alpha, n_sim, seed = seed)
  setDT(out)
  C1 <- out$C1; C2 <- out$C2
  raw[[nm]] <<- data.table(model = nm, C1 = C1, C2 = C2)
  took <- (C1 + C2) > takeoff_frac * N_tot          # major-outbreak realisations
  data.table(
    model      = nm,
    p_takeoff  = mean(took),
    mean_C1    = round(mean(C1), 1),
    mean_C2    = round(mean(C2), 1),
    ratio_mean = mean(C2) / mean(C1),
    mode_C1    = if (any(took)) round(median(C1[took]), 1) else NA_real_,
    mode_C2    = if (any(took)) round(median(C2[took]), 1) else NA_real_,
    ratio_mode = if (any(took)) median(C2[took]) / median(C1[took]) else NA_real_)
}))
summ[, ratio_gap := round(ratio_mode - ratio_mean, 3)]
summ[, ratio_mean := round(ratio_mean, 3)]
summ[, ratio_mode := round(ratio_mode, 3)]
summ[, bimodal := p_takeoff > 0.02 & p_takeoff < 0.98]
summ[, p_takeoff := round(p_takeoff, 3)]

cat(sprintf("\n=== C2/C1 ratio at mean vs mode  (beta=%.2f, alpha=%.2f, n_sim=%d) ===\n",
            beta, alpha, n_sim))
print(summ)
cat("\nbimodal = TRUE  => mode != mean => kernel and mean fits diverge => VE moves.\n")

dir.create("output", showWarnings = FALSE)
fwrite(summ, "output/mode_vs_mean_ratio.csv")

# Optional: final-size (C1) histograms per model, to see the bimodality.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  d <- rbindlist(raw)
  d[, model := factor(model, levels = names(models))]
  p <- ggplot(d, aes(C1)) +
    geom_histogram(bins = 50, fill = "grey40") +
    facet_wrap(~ model, scales = "free") +
    theme_bw(base_size = 12) +
    labs(x = "unvaccinated final size C1", y = "count",
         title = sprintf("Final-size distribution per model (beta=%.2f, alpha=%.2f)",
                         beta, alpha))
  ggsave("output/mode_vs_mean_hist.png", p, width = 11, height = 6, dpi = 130)
  cat("Wrote output/mode_vs_mean_ratio.csv and output/mode_vs_mean_hist.png\n")
}

# ---------------------------------------------------------------------------
# Does the vaccine ratio (C2/C1) depend on OUTBREAK SIZE?  A systematic
# gradient (ratio -> 1 as outbreaks shrink, because small outbreaks are
# degree/frailty-biased toward high-exposure individuals the vaccine can't
# protect) would MOVE the VE distribution. A flat cloud => small outbreaks
# only add noise (widen, not shift). Pooled ratio per size bin is robust to
# the 0/0 problem in tiny outbreaks.
# ---------------------------------------------------------------------------
size_ratio <- function(nm, nbin = 12) {
  d <- copy(raw[[nm]]); d[, total := C1 + C2]
  d <- d[total > 0]
  br <- unique(quantile(d$total, seq(0, 1, length.out = nbin + 1)))
  if (length(br) < 3) return(data.table())          # too few distinct sizes
  d[, bin := cut(total, breaks = br, include.lowest = TRUE)]
  d[, .(n = .N, med_size = median(total),
        pooled_ratio = sum(C2) / sum(C1)), by = bin][order(med_size)]
}

for (nm in intersect(c("network", "sir_sus_frailty"), names(raw))) {
  dn <- copy(raw[[nm]]); dn[, total := C1 + C2]
  rho <- suppressWarnings(cor(dn$total, ifelse(dn$C1 > 0, dn$C2 / dn$C1, NA),
                              method = "spearman", use = "complete.obs"))
  tab <- size_ratio(nm)
  cat(sprintf("\n=== %s: pooled C2/C1 by outbreak-size bin  (Spearman size~ratio = %.3f) ===\n",
              nm, rho))
  if (nrow(tab)) print(tab[, .(med_size = round(med_size, 1), n,
                               pooled_ratio = round(pooled_ratio, 3))])
  cat("  (negative Spearman + rising pooled_ratio at small sizes => real mechanism, VE shifts)\n")
}

if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(raw[["network"]])) {
  dn <- copy(raw[["network"]]); dn[, total := C1 + C2]
  dn <- dn[C1 > 0][, ratio := C2 / C1]
  trend <- size_ratio("network")
  ps <- ggplot(dn, aes(total, ratio)) +
    geom_point(alpha = 0.12, size = 0.6) +
    { if (nrow(trend)) geom_line(data = trend, aes(med_size, pooled_ratio),
                                 colour = "firebrick", linewidth = 1) } +
    geom_hline(yintercept = alpha, linetype = "dashed", colour = "steelblue") +
    geom_hline(yintercept = 1,     linetype = "dotted", colour = "grey40") +
    coord_cartesian(ylim = c(0, 1.6)) +
    theme_bw(base_size = 12) +
    labs(x = "outbreak size (C1 + C2)", y = "ratio  C2 / C1",
         title = "Network: does the vaccine ratio depend on outbreak size?",
         subtitle = "red = pooled ratio per size bin;  dashed = alpha,  dotted = 1 (no effect)")
  ggsave("output/network_size_vs_ratio.png", ps, width = 8, height = 5, dpi = 130)
  cat("Wrote output/network_size_vs_ratio.png\n")
}

# ---------------------------------------------------------------------------
# sir_split_effect: two non-mixing compartments with DIFFERENT vaccine effect
# (vaccinated susceptibility = alpha in A, 0.5/alpha in B), sizes split
# frac_A / (1-frac_A). One initial infection in a random person => it lands in
# compartment A w.p. frac_A, else B; only the seeded compartment ignites.
# Stratify-and-pool: run each seed-compartment, pool by compartment
# probability (exact mixture, variance-reduced vs a genuine random seed).
#
# Here mode != mean in RATIO (not just level): the MODE is the typical
# A-outbreak (its ratio), but the MEAN blends A and B outbreaks, which have
# different ratios. So a kernel(mode) vs mean fit would move VE via the ratio.
# ---------------------------------------------------------------------------
split_beta  <- 2.0     # reliably supercritical so each seeded compartment takes off
split_alpha <- 0.35    # => vac sus 0.35 in A (protective), 0.5/0.35=1.43 in B (harmful)
frac_A      <- 0.75
f_split     <- 0.5     # vaccinated fraction within each compartment

run_split_effect <- function(beta, alpha, n_sim, seed = 1) {
  N   <- N_tot
  N_A <- round(frac_A * N); N_B <- N - N_A
  Ngrp <- c(round((1 - f_split) * N_A), N_A - round((1 - f_split) * N_A),   # unvacA, vacA
            round((1 - f_split) * N_B), N_B - round((1 - f_split) * N_B))   # unvacB, vacB
  sus  <- c(1, alpha, 1, 0.5 / alpha)
  mm <- matrix(0, 4, 4)
  mm[1:2, 1:2] <- N / N_A          # within-compartment weight = sum(N)/N_comp
  mm[3:4, 3:4] <- N / N_B          # so each keeps R0 = beta/gamma
  tp <- seq(1, base$t_star, 1)
  run_seed <- function(gseed, nsim, sd, comp) {
    I0 <- integer(4); I0[gseed] <- 1L
    out <- run_stoch_cd_dust(mm, beta = beta, N = Ngrp, t = base$t_star, I_ini = I0,
                             susceptibility = sus, gamma = base$gamma, dt = base$dt,
                             timepoints = tp, n_sim = nsim, cores = cores, seed = sd)
    setDT(out); fin <- out[time == base$t_star]
    data.table(comp = comp, N_comp = if (comp == "A") N_A else N_B,
               C1 = fin$C1 + fin$C3, C2 = fin$C2 + fin$C4)   # pooled arms across comps
  }
  nA <- round(frac_A * n_sim)
  rbindlist(list(run_seed(1L, nA,          seed,      "A"),   # seed in unvac of A
                 run_seed(3L, n_sim - nA,  seed + 1L, "B")))  # seed in unvac of B
}

sp <- run_split_effect(split_beta, split_alpha, n_sim)
sp[, total := C1 + C2]
sp[, took  := total > 0.2 * N_comp]                # took off within its compartment

# Per-compartment ratio among takeoffs (case-pooled, robust to 0/0).
comp_tab <- sp[(took), .(p = .N / n_sim, mean_C1 = round(mean(C1), 1),
                         mean_C2 = round(mean(C2), 1),
                         ratio = round(sum(C2) / sum(C1), 3)), by = comp][order(comp)]
ratio_mode <- comp_tab[comp == "A", ratio]                       # dominant outcome = A
ratio_mean <- sp[(took), sum(C2) / sum(C1)]                      # blend across comps

cat(sprintf("\n=== sir_split_effect  (beta=%.2f, alpha=%.2f -> vac sus A=%.2f, B=%.2f; %d/%d split) ===\n",
            split_beta, split_alpha, split_alpha, 0.5 / split_alpha,
            round(100 * frac_A), round(100 * (1 - frac_A))))
print(comp_tab)
cat(sprintf("ratio at MODE (typical A outbreak) = %.3f  ->  VE_mode = %.3f\n",
            ratio_mode, 1 - ratio_mode))
cat(sprintf("ratio at MEAN (A+B blend)          = %.3f  ->  VE_mean = %.3f\n",
            ratio_mean, 1 - ratio_mean))
cat(sprintf("=> ratio gap (mean - mode) = %.3f: the compartment mixture MOVES VE via the ratio\n",
            ratio_mean - ratio_mode))
fwrite(comp_tab, "output/split_effect_ratio.csv")

# ---------------------------------------------------------------------------
# Model x functional variance table. For each model, how much does one
# realisation tell us about the next? Compare the BETWEEN-realisation variance
# of a functional to the i.i.d. (independent-units) benchmark -- the variance
# it WOULD have if the arms were binomial samples of the same mean attack rate.
# That ratio is the design effect = N / N_eff:
#   ~1  => interference adds nothing; one sample is as good as N indep units.
#   >>1 => shared epidemic inflates the variance; effective N collapses.
# Two functionals: the arm RATIO (log risk ratio; the shared force of infection
# CANCELS) and the arm DIFFERENCE (risk difference; it does NOT). Conditioned
# on takeoff so both are well-defined (identification given an outbreak).
# ---------------------------------------------------------------------------
functional_de <- function(dt, N1, N2, takeoff_frac = 0.05) {
  d <- copy(dt); d[, tot := C1 + C2]
  d <- d[tot > takeoff_frac * (N1 + N2) & C1 > 0 & C2 > 0]
  if (nrow(d) < 20) return(NULL)
  p1 <- mean(d$C1) / N1; p2 <- mean(d$C2) / N2          # mean attack rates
  RD  <- d$C1 / N1 - d$C2 / N2                          # risk difference
  lRR <- log((d$C2 / N2) / (d$C1 / N1))                 # log risk ratio (vac/unvac)
  Viid_RD  <- p1 * (1 - p1) / N1 + p2 * (1 - p2) / N2   # binomial benchmarks
  Viid_lRR <- (1 - p1) / (N1 * p1) + (1 - p2) / (N2 * p2)
  de_ratio <- var(lRR) / Viid_lRR
  de_diff  <- var(RD)  / Viid_RD
  data.table(n_take = nrow(d), p1 = round(p1, 3), p2 = round(p2, 3),
             DE_ratio = round(de_ratio, 1), DE_diff = round(de_diff, 1),
             diff_over_ratio = round(de_diff / de_ratio, 1),
             Neff_ratio = round((N1 + N2) / de_ratio),
             Neff_diff  = round((N1 + N2) / de_diff))
}

N1 <- base$N_cont; N2 <- base$N_vac
func_dt <- rbindlist(lapply(names(raw), function(nm)
  cbind(model = nm, functional_de(raw[[nm]], N1, N2))), fill = TRUE)
sp_pooled <- functional_de(sp[, .(C1, C2)], N1, N2)     # split-effect mixture
if (!is.null(sp_pooled)) func_dt <- rbind(func_dt, cbind(model = "sir_split_effect", sp_pooled))

cat(sprintf("\n=== Model x functional design effect (N=%d; DE = between-var / i.i.d.-var = N/N_eff) ===\n",
            N1 + N2))
print(func_dt[, .(model, DE_ratio, DE_diff, diff_over_ratio, Neff_ratio, Neff_diff)])
cat("  DE ~ 1: functional robust to interference; DE >> 1: variance inflated, N_eff collapses.\n")
cat("  Standard models: ratio (shared FOI cancels) well identified, difference not.\n")
cat("  sir_split_effect: even the ratio is inflated (we broke the cancellation).\n")
fwrite(func_dt, "output/model_functional_variance.csv")
