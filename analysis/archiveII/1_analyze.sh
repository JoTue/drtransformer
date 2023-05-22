#!/bin/bash

#SBATCH --job-name=1_a
#SBATCH --partition=skylake_0096
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH --output=log/1_%x_%j.out
#SBATCH --error=log/1_%x_%j.err

NTHREADS=1

DATASET=$1  # e.g. "tRNA.db"
ID1=$2
ID2=$3
IDS_PER_JOB=$4

LF=log

if [ $ID1 -eq 0 ] && [ $ID2 -eq 0 ] # analyze complete dataset
then
    ID1=1
    ID2=$(grep -c '^>' data/$DATASET)
fi

# analyze ID1 to ID2 (both inclusive)
seq $ID1 $IDS_PER_JOB $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "bin/analyze.py data/$DATASET --id1 {} --id-end $ID2 --ids-per-job $IDS_PER_JOB --tmp analysis_tmp"