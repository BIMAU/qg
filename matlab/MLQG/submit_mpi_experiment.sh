#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --ntasks=4
#SBATCH --nodes=4
#SBATCH -p short

module load MCR/R2018a
export MCR_CACHE_ROOT=`mktemp -d /scratch-local/mcr.XXXXXX`

cd /home/emulder/Projects/qg/matlab/MLQG

mkdir -p data/experiments/logdir
rm data/experiments/logdir/*

date=`date +%m%d%y-%H%M`
procs=4  # Actually the number of nodes. Matlab does threading on a
         # node.

echo "MLQG experiment"
echo " #procs = " $procs

srun -n $procs ./mpi_experiment
