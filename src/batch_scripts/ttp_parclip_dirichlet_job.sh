#!/bin/sh
#SBATCH -N 3
#SBATCH --ntasks-per-node 5
#SBATCH --output ttp_parclip_dirichlet.out
#SBATCH --qos himem
mpiexec -np 5 lb run_ttp_parclip_dirichlet.sh
