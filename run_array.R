# Slurm array runner: each task picks a chunk of `configs` and runs
# `run_from_config` over it with mclapply, writes results to RDS keyed by
# the array task id.
#
# Submit (example, adjust partition/time/cpus):
#   sbatch --array=1-20 --cpus-per-task=10 run_array.sh
# where run_array.sh wraps `Rscript run_array.R`.

library(dplyr)
library(data.table)
library(glue)

source("det_model.R")
source("stoch_model.R")
source("utils.R")
source("estimate_from_data.R")  # builders + run_from_config + esteimate_VE_from_fs

RhpcBLASctl::omp_set_num_threads(1)
setDTthreads(1)

# ---------------------------------------------------------------------------
# Build the config grid
#   1 linear model
#   1 sir model
#   10 frailty configs (vary allocation_seed)
#   10 contact matrices x 10 allocations = 100 network configs
# Total: 112 configs.
# ---------------------------------------------------------------------------

N_cont <- 100; N_vac <- 100; X_cont <- 50; X_vac <- 25; t <- 8

base <- list(
    N        = c(N_cont, N_vac),
    t        = t,
    dt       = 0.01,
    timepoints = seq(1, t, 1),
    n_sim    = 1000,
    cores    = 4,                                    # dust threads per call
    # Stochastic frozen-field EATE: eate_n_vac is the number of outer
    # vac re-allocations per call, eate_n_rep is the number of dust
    # replicates per allocation.
    eate_n_vac = 10, eate_n_rep = 20,
    N_vac    = N_vac, X_vac = X_vac, N_cont = N_cont, X_cont = X_cont,
    n_cores  = 4,                                    # calc_ve_from_pars mclapply
    beta_ini = 1.5, alpha_ini = 0.5,
    n = 50, burn_in = 50, scale = c(0.1, 0.05)
)

configs <- list()

# 1 linear
configs[[length(configs)+1]] <- modifyList(base, list(
    name      = "linear",
    generator = purrr::partial(run_stoch_linear_dust,
                               N=c(N_cont, N_vac), t=t, dt=0.01,
                               timepoints=seq(1, t, 1), n_sim=1000, cores=4),
    eate_func = get_EATE_linear,
    beta_ini  = 0.05, alpha_ini = 0.5,
    n = 500, burn_in = 500, scale = c(0.02, 0.05)))

# 1 sir (homogeneous mixing). eate_func must match the nonlinear simulator:
# get_EATE (contact-dependent), not get_EATE_linear. gamma=1 is the project
# convention and must be set on both the simulator and the EATE function
# (otherwise their defaults diverge: stoch_cd_dust=1/3, get_EATE=0.2).
configs[[length(configs)+1]] <- modifyList(base, list(
    name      = "sir",
    generator = purrr::partial(run_stoch_cd_dust,
                               matrix(rep(1, 4), nrow=2),
                               I_ini=c(10,10), N=c(N_cont, N_vac),
                               gamma=1,
                               t=t, dt=0.01,
                               timepoints=seq(1, t, 1), n_sim=1000, cores=4),
    eate_func = purrr::partial(get_EATE, gamma=1)))

# 10 frailty configs — vary allocation_seed only
for (allocation_seed in 1:10) {
    configs[[length(configs)+1]] <- modifyList(base, list(
        name        = glue("frailty_a{allocation_seed}"),
        model_type  = "frailty",
        sd          = 0,
        sd_trans    = 0.3,
        n_frailty   = 10,
        I_ini       = c(10, 10),
        gamma       = 1,
        allocation_seed = allocation_seed))
}

# 10 contact matrices x 10 allocations = 100 network configs
for (network_seed in 1:10) {
    for (allocation_seed in 1:10) {
        configs[[length(configs)+1]] <- modifyList(base, list(
            name        = glue("network_n{network_seed}_a{allocation_seed}"),
            model_type  = "network",
            pl_alpha    = 3,
            mean_k      = 6,
            I_ini       = 2,
            gamma       = 1,
            # 10 outer allocation_seeds already give us 10 vac realisations
            # per network; no need for more inside get_stoch_eate_network.
            eate_n_vac      = 1,
            network_seed    = network_seed,
            allocation_seed = allocation_seed))
    }
}

message(glue("Built {length(configs)} configs."))

# ---------------------------------------------------------------------------
# Slurm array dispatch
# ---------------------------------------------------------------------------

N_nodes        <- 20      # match --array=1-N_nodes in the sbatch call
cores_per_node <- 10      # outer mclapply (each call uses `cores` inside dust)

id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (id == "") id <- 1
id <- as.numeric(id)

chunk2  <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))
set.seed(1)
chunked <- chunk2(sample(seq_along(configs)), N_nodes)
set.seed(NULL)

ids <- as.integer(chunked[[id]])
message(glue("Task {id}: running configs {paste(ids, collapse=',')}"))

run_id <- function(i) {
    message(glue("Start {i}: {configs[[i]]$name}"))
    out <- tryCatch(run_from_config(configs[[i]]),
                    error = function(e) { message("Error: ", conditionMessage(e)); NULL })
    message(glue("End {i}"))
    list(name = configs[[i]]$name, config = configs[[i]], result = out)
}

results <- parallel::mclapply(ids, run_id,
                              mc.cores = cores_per_node,
                              mc.preschedule = FALSE)

out_dir <- "output/array_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(out_dir, glue("results_{id}.RDS")))
message(glue("Wrote {out_dir}/results_{id}.RDS"))
