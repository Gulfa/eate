#!/bin/bash
#SBATCH --job-name=eate_fit_array
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --output=output/fit_array_results/slurm_%A_%a.out
#SBATCH --error=output/fit_array_results/slurm_%A_%a.err

mkdir -p output/fit_array_results

Rscript run_fit_array.R
