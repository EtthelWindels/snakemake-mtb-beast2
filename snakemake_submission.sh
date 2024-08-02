#!/bin/bash

#SBATCH -n 1
#SBATCH --time=360:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --tmp=1000
#SBATCH --job-name=snakemake
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err


# Note: snakemake-phylo-beast2/slurm should be in ~/.config/snakemake

source activate snakemake
module load jdk
module load stack/2024-06 r/4.3.2
cd path_to_snakemake-mtb-beast2
snakemake --profile slurm --unlock
snakemake --cores 1000 --profile slurm --use-conda --keep-going --ignore-incomplete --rerun-triggers mtime