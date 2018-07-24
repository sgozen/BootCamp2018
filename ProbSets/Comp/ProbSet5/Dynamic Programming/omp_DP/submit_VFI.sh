#!/bin/bash

#SBATCH --job-name=DP
#SBATCH --output=VFI1.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sandyb
#SBATCH --account=osmlab

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the process
./VFI
