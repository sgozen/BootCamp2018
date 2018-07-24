#!/bin/bash -l

#SBATCH --job-name=pi
#SBATCH --output=pi.out
#SBATCH --error=pi.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --partition=sandyb
#SBATCH --time=01:00:00

module load openmpi

### MPI executable
mpiexec -np $SLURM_NTASKS ./pi.exe
