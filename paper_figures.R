# Paper figures, in paper order, apart from the fit-to-data figures
# (run_fit_array.R / analyse_fit_array.R, which need the slurm array).
#
#   Rscript paper_figures.R
#
# Each figure has its own script, which writes its working outputs to
# output/; the panel used in the paper is then copied to output/paper/ under
# the same file name.
# run_paper.R is the older, exploratory pipeline and is not needed here.

dir.create("output/paper", recursive = TRUE, showWarnings = FALSE)

paper_figure <- function(n, script, file) {
  message(sprintf("Figure %d: %s", n, script))
  source(script, local = new.env())
  stopifnot(file.copy(file, file.path("output/paper", basename(file)), overwrite = TRUE))
}

# Fig 1: HR and CRR/CIR under mean-1 gamma frailty, linear (exposure-only) model
paper_figure(1, "linear_frailty_hr_crr.R", "output/attack_rate_hetero_gamma.png")

# Fig 2: attack rate and VE over time, SIR with a fixed leaky vaccine
paper_figure(2, "sir_attack_rate_ve.R", "output/attack_rate_nl.png")

# Fig 3: VE(t) by alpha across linear, SIR, frailty and network models
# (slow: sized for the server)
paper_figure(3, "ve_by_alpha.R", "output/ve_by_alpha.png")
