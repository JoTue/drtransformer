#!/bin/bash

#SBATCH --job-name=0_d
#SBATCH --partition=skylake_0096
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --output=log/0_%x_%j.out
#SBATCH --error=log/0_%x_%j.err

NTHREADS=16

DATASET=$1  # e.g. "tRNA.db"
ID1=$2
ID2=$3
IDS_PER_JOB=$4

LF=log

if [ $ID1 -eq 0 ] && [ $ID2 -eq 0 ] # compute complete dataset
then
    ID1=1
    ID2=$(grep -c '^>' data/$DATASET)
fi

# compute ID1 to ID2 (both inclusive)
seq $ID1 $IDS_PER_JOB $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "sbatch -o $LF/d_%A.out -e $LF/d_%A.err -J d-{}-$DATASET -p skylake_0096 --nice=0 bin/drtransformer.py data/$DATASET --id1 {} --id-end $ID2 --ids-per-job $IDS_PER_JOB"
