#!/bin/sh
#SBATCH --output=mpi.out
#SBATCH --ntasks=4
#SBATCH --partition=cpar
perf stat -e instructions,cycles mpirun -np 4 bin/k_means 10000000 4
