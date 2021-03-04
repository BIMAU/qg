#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: ./compile_mfile main.m build_dir"
    exit
fi

mkdir -p $2

module load MATLAB/2018a
module load MCR/R2018a

export MCR_CACHE_ROOT=`mktemp -d /scratch-local/mcr.XXXXXX`

mcc -R -singleCompThread -v -C -m \
    -a ~/local/matlab \
    -a ~/Projects/ESN/matlab/ESN.m \
    -d $2 \
    $1

#mcc -v -C -m \
#    -a ~/local/matlab \
#    -a ~/Projects/ESN/matlab/ESN.m \
#    $1


cp -v mpi_experiment.cc $2/.
cp -v submit_mpi_experiment.sh $2/.

cd $2

mpicxx -g -Wall mpi_experiment.cc -o mpi_experiment -lmpi

sbatch submit_mpi_experiment.sh

cd -
