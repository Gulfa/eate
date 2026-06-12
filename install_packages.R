# Install all dependencies for the eate project.
#   Rscript install_packages.R

options(repos = c(CRAN = "https://cran.rstudio.com"))

cran_pkgs <- c(
  # data wrangling / plotting
  "dplyr", "data.table", "tidyr", "purrr", "glue", "zoo",
  "ggplot2", "ggtext", "cowplot", "scales", "latex2exp",
  # numerics / stats
  "adaptivetau", "deSolve", "Pareto", "boot", "distcrete",
  "mcmc", "mvtnorm", "coda",
  # parallel / system
  "RhpcBLASctl", "remotes",
  # legacy odin (CRAN). odin2 + dust2 below via GitHub.
  "odin",
  # tests
  "testthat"
)

to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
} else {
  message("All CRAN packages already installed.")
}

# odin2 + dust2 are only on GitHub. Memory note: prefer odin2 + dust2 with
# n_particles over mclapply for stochastic replication in this project.
gh_pkgs <- list(
  dust2     = "mrc-ide/dust2",
  odin2     = "mrc-ide/odin2",
  odin.dust = "mrc-ide/odin.dust"
)
for (pkg in names(gh_pkgs)) {
  if (!pkg %in% rownames(installed.packages())) {
    message("Installing ", pkg, " from ", gh_pkgs[[pkg]])
    remotes::install_github(gh_pkgs[[pkg]])
  } else {
    message(pkg, " already installed.")
  }
}

message("Done.")
