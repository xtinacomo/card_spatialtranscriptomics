#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 96:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#Load required modules
module purge
module load apptainer
module load snakemake/7.7.0

# Optional: Only clone profile if not already cloned
if [ ! -d "snakemake_profile" ]; then
  git clone https://github.com/NIH-HPC/snakemake_profile.git
fi

# Pull the containers
mkdir -p envs
#apptainer pull --disable-cache envs/snapatac2.sif oras://quay.io/adamcatchingdti/snapatac2
#apptainer pull --disable-cache envs/single_cell_gpu.sif oras://quay.io/adamcatchingdti/single_cell_gpu:0.8
#apptainer pull --disable-cache envs/decoupler.sif oras://quay.io/adamcatchingdti/decoupler.sif:0.9

#apptainer pull --disable-cache envs/scenicplus.sif docker://litd/docker-scenicplus:latest 

# Load singularity
module load singularity

# Bind external directories on Biowulf
#. /usr/local/current/singularity/app_conf/sing_binds

# Make sure script folder permissions are correct
#chmod +x scripts/*.sh

# Create logs directory if it doesn't exist
mkdir -p logs

# Run Snakemake with Singularity support and helpful flags
snakemake \
  --snakefile Snakefile \
  --profile snakemake_profile \
  --configfile config.yaml \
  --use-singularity \
  --cores all \
  --printshellcmds \
  --reason \
  --latency-wait 60 \
  --rerun-incomplete

