#!/bin/sh
#SBATCH --output rf00037_test_shuffle.out
#SBATCH --qos himem
#SBATCH --mem=250G
#srun lb rf00037_run_shuffle.sh
srun -n 3 lb rf00037_run_shuffle_seqbased.sh
