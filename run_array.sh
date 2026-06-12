#!/bin/bash
#SBATCH --job-name=eate_array
#SBATCH --array=1-20
#SBATCH --cpus-per-task=10
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --output=output/array_results/slurm_%A_%a.out
#SBATCH --error=output/array_results/slurm_%A_%a.err

mkdir -p output/array_results

Rscript run_array.R
