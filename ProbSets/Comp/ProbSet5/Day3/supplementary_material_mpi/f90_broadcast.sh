#!/bin/bash -l

#SBATCH --ntasks=16

#SBATCH --time=00:02:00

#SBATCH --output=f90.broadcast.out
#SBATCH --error=f90.broadcast.err


### MPI executable
mpirun -np $SLURM_NTASKS ./broadcast.exec
