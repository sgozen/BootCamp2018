#!/bin/bash -l

#SBATCH --ntasks=8


#SBATCH --time=00:02:00

#SBATCH --output=f90.est_pi8.out
#SBATCH --error=f90.est_pi8.err


### MPI executable
mpirun -np $SLURM_NTASKS ./est_pi.exec
