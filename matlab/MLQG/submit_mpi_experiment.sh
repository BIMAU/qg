#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=8
#SBATCH --tasks-per-node=8
#SBATCH -p normal

module load MCR/R2018a
export MCR_CACHE_ROOT=`mktemp -d /scratch-local/mcr.XXXXXX`

cd /home/emulder/Projects/qg/matlab/MLQG

# mkdir -p data/experiments/logdir
# rm data/experiments/logdir/*

srun -n 8 ./mpi_experiment
