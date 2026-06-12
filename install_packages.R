# Install all dependencies for the eate project.
#   Rscript install_packages.R

options(repos = c(CRAN = "https://cran.rstudio.com"))

# dust2 / odin2 need C++17 (if constexpr, std::is_invocable, structured
# bindings). Some HPC R builds (e.g. R 4.2 foss-2022a) leave CXX17 unset
# in Makeconf, so package compilation falls back to gnu++14 and fails.
# We persist the override to ~/.R/Makevars because R_MAKEVARS_USER does
# not always propagate to install subprocesses on HPC. The -std= flag is
# folded into CXX17 itself so it survives even when CXX17STD is ignored.
ensure_cxx17 <- function() {
  lines <- c(
    "CXX17 = g++ -std=gnu++17",
    "CXX17STD = -std=gnu++17",
    "CXX17FLAGS = -O2 -fPIC -fopenmp",
    "CXX17PICFLAGS = -fPIC"
  )
  r_dir <- path.expand("~/.R")
  if (!dir.exists(r_dir)) dir.create(r_dir, recursive = TRUE)
  mv <- file.path(r_dir, "Makevars")
  current <- if (file.exists(mv)) readLines(mv) else character()
  to_add <- setdiff(lines, current)
  # Drop any prior conflicting CXX17 lines so the new values win.
  if (length(to_add) > 0) {
    pruned <- current[!grepl("^\\s*CXX17(STD|FLAGS|PICFLAGS)?\\s*=", current)]
    writeLines(c(pruned, lines), mv)
    message("Updated ", mv, " with CXX17 settings")
  } else {
    message(mv, " already has CXX17 settings")
  }
  # Belt-and-suspenders: also point R_MAKEVARS_USER at the same file for
  # subprocesses that might honour it.
  Sys.setenv(R_MAKEVARS_USER = mv)
  invisible(mv)
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
