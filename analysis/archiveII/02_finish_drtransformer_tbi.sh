#!/bin/bash

#SBATCH --job-name=02_finish
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --output=log/02_f.out
#SBATCH --error=log/02_f.err

NTHREADS=16

DATASET=$1  # e.g. "tRNA.db"
DATASET_name="$(basename $DATASET .db)"
ID1=$2
ID2=$3

LF=log

if [ $ID1 -eq 0 ] && [ $ID2 -eq 0 ] # compute complete dataset
then
    ID1=1
    ID2=10000
fi

# compute ID1 to ID2 (both inclusive)
bin/get_missing_ids_from_file.py $DATASET --id1 $ID1 --id2 $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "bin/drtransformer.py data/$DATASET --id {}"
