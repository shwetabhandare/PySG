#!/bin/sh
#SBATCH --output lin28_test_shuffle.out
#SBATCH --qos himem
#SBATCH --mem=250G
#SBATCH --time 48:00:00
#srun lb lin28_run_shuffle.sh
srun -n 3 lb lin28_run_shuffle_seqbased.sh
