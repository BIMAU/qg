#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: ./compile_mfile main.m"
fi

module load matlab
module load mcr

export MCR_CACHE_ROOT=`mktemp -d /scratch-local/mcr.XXXXXX`

mcc -R -singleCompThread -v -C -m $1
