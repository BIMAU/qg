#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: ./compile_mfile main.m"
fi

module load MATLAB/2018a
module load MCR/R2018a

export MCR_CACHE_ROOT=`mktemp -d /scratch-local/mcr.XXXXXX`

mcc -R -singleCompThread -v -C -m \
    -a ~/local/matlab \
    -a ~/Projects/ESN/matlab/ESN.m \
    $1
