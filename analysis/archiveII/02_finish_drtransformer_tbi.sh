#!/bin/bash

#SBATCH --job-name=02_f
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --output=log/02_%x_%j.out
#SBATCH --error=log/02_%x_%j.err

NTHREADS=16

DATASET=$1  # e.g. "tRNA.db"
DATASET_name="$(basename $DATASET .db)"
ID1=$2
ID2=$3

# compute ID1 to ID2 (both inclusive)
bin/get_missing_ids_from_file.py $DATASET --id1 $ID1 --id2 $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "bin/drtransformer.py data/$DATASET --id {}"
