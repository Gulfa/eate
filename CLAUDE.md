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
- `tests/testthat/test-net-sir-events.R` — the event-driven network simulator and the `run_stoch_network()` engine switch (needs `cpp11`; compiles `net_sir_events.cpp` on first use).
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

## Model compile cache (`odin_cache.R`)

The six odin2/dust2 models take ~3 minutes to compile from scratch. `odin_cache.R`
points dust2 at a persistent build directory (`.odin_cache/`, gitignored) via
`DUST_WORKDIR_ROOT`, so a fresh session reuses the compiled shared libraries and
starts in ~6s instead. It is loaded automatically by the project `.Rprofile` and
by `model.R` / `det_model.R` / `stoch_model.R`.

Nothing needs to be invalidated by hand: the cache key is a hash of the generated
C++, so editing a model file recompiles it and reverting the edit gets the old
build back for free. Compiles are locked per model, so the slurm arrays can start
cold without racing each other.

```r
odin_cache_warm()          # precompile everything (do this before submitting an array)
odin_cache_status()        # what is cached, when it was built, how big
odin_cache_clean(days=30)  # reclaim disk; 0 clears the cache
```

Run with `R --vanilla` or set `EATE_ODIN_CACHE_DISABLE=1` to bypass it. The two
odin v1 models (`det_mod_cd.R`, `det_mod_ncd.R`) are not cached — odin v1 rewrites
its generated C every time — but they are only ~5s together.

## Network simulator engine (`net_sir_events.cpp`)

`run_stoch_network()` has two interchangeable backends, selected by `engine=`:

| engine | implementation | notes |
|---|---|---|
| `"events"` *(default)* | `net_sir_events.cpp` via cpp11 | exact event-driven (first-passage percolation), no `dt` bias, O(E log V), OpenMP over realisations |
| `"dust"` | `stoch_mod_adj.R` via `run_stoch_adj()` | discrete-time tau-leaping; what every result before this switch used |

Both simulate the same continuous-time process and return the same
`time, sim, C1, C2` columns on the same `beta` scale. The event engine ignores
`dt` (it has no time step) and is ~50–250× faster at the `dt` the fit needs.

Switch globally, without editing code:

```sh
EATE_NETWORK_ENGINE=dust Rscript run_fit_array.R    # env var wins
```
```r
options(eate.network_engine = "dust")               # or per session
run_stoch_network(..., engine = "dust")             # or per call
```

`run_fit_array.R` resolves the engine once into `net_engine`, records it on
every network result as `network_engine`, and prebuilds the CSR graph
(`cfg$.csr`, via `adj_to_csr()`) once per config alongside `cfg$.adj`.

Two things the event engine does **not** cover:

- **Per-node trajectories.** It returns infection times only, no S/I/R paths,
  so `get_stoch_eate_network()` (which integrates a per-node force of
  infection) always uses dust. Only the fit/simulator side switches.
- **Time-varying rates.** Validity requires scalar `beta`, constant `gamma`,
  and `waning = 0`.

Index cases are excluded from `C1`/`C2` by both engines (`stoch_mod_adj.R` has
`initial(C[]) <- 0`); `run_stoch_network_events(count_seeds = TRUE)` includes
them if you want incidence with the seeds.

`validate_net_sir_events.R` compares the two across a `dt` ladder — the dust
bias must shrink towards the event engine as `dt → 0`.

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
