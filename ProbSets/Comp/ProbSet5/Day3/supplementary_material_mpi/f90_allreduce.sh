#!/bin/bash -l
#SBATCH --ntasks=10

#SBATCH --time=00:02:00

#SBATCH --output=f90.allreduce.out
#SBATCH --error=f90.allreduce.err


### MPI executable
mpirun -np $SLURM_NTASKS ./allreduce.exec
