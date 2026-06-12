# Install all dependencies for the eate project.
#   Rscript install_packages.R

options(repos = c(CRAN = "https://cran.rstudio.com"))

# dust2 / odin2 need C++17 (if constexpr, std::is_invocable, structured
# bindings). Some HPC R builds (e.g. R 4.2 foss-2022a) leave CXX17 unset
# in Makeconf, so package compilation falls back to gnu++14 and fails.
# Force C++17 by pointing R_MAKEVARS_USER at a temp Makevars for this run.
# Persist by appending the same three lines to ~/.R/Makevars if desired.
ensure_cxx17 <- function() {
  mv <- tempfile("Makevars_cxx17_")
  writeLines(c(
    "CXX17 = g++",
    "CXX17STD = -std=gnu++17",
    "CXX17FLAGS = -O2 -fPIC -fopenmp"
  ), mv)
  old <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
  Sys.setenv(R_MAKEVARS_USER = mv)
  message("Forcing C++17 via temp Makevars: ", mv)
  invisible(old)
}
ensure_cxx17()

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
