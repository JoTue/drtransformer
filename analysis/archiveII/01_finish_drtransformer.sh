#!/bin/bash

#SBATCH --job-name=01_finish
#SBATCH --partition=skylake_0096
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --output=log/01_f.out
#SBATCH --error=log/01_f.err

NTHREADS=16

DATASET=$1  # e.g. "tRNA.db"
DATASET_name="$(basename $DATASET .db)"
ID1=$2
ID2=$3
print_only_id=${4:-0}

LF=log

if [ $ID1 -eq 0 ] && [ $ID2 -eq 0 ] # compute complete dataset
then
    ID1=1
    ID2=$(grep -c '^>' data/$DATASET)
fi

# compute ID1 to ID2 (both inclusive)
# seq $ID1 $IDS_PER_JOB $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "bin/drtransformer.py data/$DATASET --id1 {} --id-end $ID2 --ids-per-job $IDS_PER_JOB"
bin/get_missing_ids.py data/$DATASET drtr_results/$DATASET_name --id1 $ID1 --id2 $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "bin/drtransformer.py data/$DATASET --id {} --print-only-id $print_only_id"
