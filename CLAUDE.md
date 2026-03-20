# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based epidemiological research project studying **Expected Average Treatment Effect (EATE)** — a causal inference framework for vaccine effectiveness (VE) estimation in contact networks. It uses stochastic simulation and counterfactual analysis to decompose direct and indirect vaccination effects.

## Tests

Tests use `testthat`. Run from the project root:

```r
# All tests
Rscript tests/run_tests.R

# Or interactively
library(testthat)
test_dir("tests/testthat")

# Single file
testthat::test_file("tests/testthat/test-stacking.R")
```

- `tests/testthat/test-model.R` — tests for `get_frailty`, `cij_NGM`, `get_beta`, `run_det_cd`, `run_det_ncd`. Requires `odin` and the odin model files; must be run from the project root (`chdir=TRUE` is set automatically).
- `tests/testthat/test-stacking.R` — tests for `regularise`, `generate_contacts`, `generate_linear_event_times`, `run_events_linear`, `run_stacking`, plus a validation test comparing the custom linear simulator against `adaptivetau`. Does **not** require odin.

**Note:** `stacking.R` contains several informal test/validation functions (`test()`, `test_nl()`, `test_linear_implementation()`, `test_sir_implementation()`) that print comparison output but have no assertions. `def_run_test_linear()` is currently broken — it references `res$stacked` which was commented out of `compare_cf`.

## Running the Code

This is a research codebase with no formal build system. All scripts are run interactively in R or RStudio:

```r
source("model.R")      # Load core simulation functions
source("stacking.R")   # Load stochastic counterfactual functions
source("run_paper.R")  # Execute full paper analysis pipeline and generate figures
```

**Required R packages:**
```r
install.packages(c("dplyr", "data.table", "ggplot2", "adaptivetau", "odin",
                   "glue", "tidyr", "cowplot", "Pareto", "scales", "zoo"))
```

The `odin` package may need to be installed from GitHub: `remotes::install_github("mrc-ide/odin")`

## Architecture

The code is organized in three layers:

### Layer 1: ODE Model Definitions (`det_mod_cd.R`, `det_mod_ncd.R`)
Written in `odin` DSL format — these are compiled at runtime via `odin::odin()`, not executed directly. They define compartmental SIR models:
- `det_mod_cd.R`: Contact-dependent model with mixing matrix C_ij, time-varying beta, and per-group susceptibility modifiers
- `det_mod_ncd.R`: Simpler non-contact-dependent variant (S→R with susceptibility modifiers)

### Layer 2: Core Functions (`model.R`)
Loads and compiles the odin models, then wraps them with higher-level R functions:
- `run_det_cd()` / `run_det_ncd()` — run the ODE models and return effect measures
- `run_frailty_cd()` / `run_frailty()` — add beta-distributed susceptibility heterogeneity
- `get_eate_frailty()` / `get_eate_network()` — causal EATE estimation
- `get_conact_matrix_pl()` — generates Pareto (power-law) contact networks
- `run_mean_field()` — mean-field approximation for stratified populations
- `get_beta()` — derives transmission rate β from basic reproduction number R₀
- `cij_NGM()` — constructs next-generation matrix

### Layer 3: Analysis (`stacking.R`, `run_paper.R`)
- `stacking.R` implements coupled/stacked tau-leaping stochastic simulations using `adaptivetau`. The "stacking" approach runs factual and counterfactual scenarios with shared random draws for variance reduction. Key functions: `run_stacking()`, `compare_cf()`, `calc_EATEs()`, `DE_sir()`, `TO_sir()`
- `run_paper.R` is the top-level analysis script that sources both files above, runs all analyses, and generates figures saved to `output/` and `article/figures/`

## Key Concepts

**Effect measures computed:**
- `VE` — Vaccine Effectiveness (1 − Risk Ratio)
- `HRR` — Hazard Rate Ratio (instantaneous)
- `CIR` — Cumulative Incidence Ratio (attack rate)
- `CRR` — Cumulative Rate Ratio
- `EATE` — Expected Average Treatment Effect (causal)

**Causal decomposition:**
- `DE` (Direct Effect) — effect of vaccinating one person on themselves
- `TO` (Total Outcome) — combined individual + indirect network effect
- EATE is estimated via leave-one-out perturbation: iteratively removing individuals and averaging causal effects

**Heterogeneity modeling:**
- Frailty: multiplicative beta-distributed susceptibility modifier (parameter `sd` controls spread)
- Stratification: population split into groups with different contact rates and susceptibilities (the `alpha` parameter)
- Network: Pareto-distributed degree distribution for realistic contact structure
