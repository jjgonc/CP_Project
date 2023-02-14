#!/bin/bash
#SBATCH --output=mpi.out
#SBATCH --ntasks=16
#SBATCH --partition=cpar
time mpirun -np 16 bin/k_means 10000000 4