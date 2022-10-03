#!/bin/bash

#SBATCH --job-name=ar_drtr
#SBATCH --partition=skylake_0096
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1

NTHREADS=16

DATASET=$1  # e.g. "tRNA.db"
ID1=$2
ID2=$3

LF=log

if [ $ID1 -eq 0 ] && [ $ID2 -eq 0 ] # compute complete dataset
then
    ID1=1
    ID2=$(grep -c '^>' data/$DATASET)
fi

# compute ID1 to ID2 (both inclusive)
seq $ID1 $IDS_PER_JOB $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "sbatch -o $LF/%A_%a.out -e $LF/%A_%a.err -J drtr-{}-$DATASET -a $ID1-$ID2 -p skylake_0096 --nice=1000 bin/drtransformer.py data/$DATASET --id1 {} --id2 $(($(($ID1 + $IDS_PER_JOB - 1)) < $ID2 ? $(($ID1 + $IDS_PER_JOB - 1)) : $ID2))"
