#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --tasks=12
#SBATCH -p short
#SBATCH --mem-per-cpu=2000

module load MCR/R2018a
export MCR_CACHE_ROOT=`mktemp -d /scratch-local/mcr.XXXXXX`

cd /home/emulder/Projects/qg/matlab/MLQG

mkdir -p data/experiments/logdir
rm data/experiments/logdir/*

date=`date +%m%d%y-%H%M`
nodes=1
procs=12

echo "MLQG experiment"
echo " #procs = " $procs

srun -n $nodes ./mpi_experiment $procs
