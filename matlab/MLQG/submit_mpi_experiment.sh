#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --tasks=20
#SBATCH --tasks-per-node=5
#SBATCH -p short

module load MCR/R2018a
export MCR_CACHE_ROOT=`mktemp -d /scratch-local/mcr.XXXXXX`

srun -n 20 ./mpi_experiment
