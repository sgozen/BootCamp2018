#!/bin/bash -l

#SBATCH --ntasks=5

#SBATCH --time=00:02:00

#SBATCH --output=f90.scatter.out
#SBATCH --error=f90.scatter.err


### MPI executable
mpirun -np $SLURM_NTASKS ./scatter.exec
