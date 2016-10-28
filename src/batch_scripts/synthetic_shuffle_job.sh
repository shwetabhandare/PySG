#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output synthetic_shuffle.out
#SBATCH --qos himem
srun lb run_synthetic_dirichlet.sh
