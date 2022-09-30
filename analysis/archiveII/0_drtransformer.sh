#!/bin/bash

#SBATCH --job-name=drtr
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --mem=20G

NTHREADS=16

DATASET=$1
ID1=$2
ID2=$3

LF=log

if [ $ID1 -eq 0 ] && [ $ID2 -eq 0 ] # compute complete dataset
then
    ID1=1
    ID2=$(grep -c '^>' data/$DATASET)
fi

# compute ID1 to ID2 (both inclusive), each as separate job
seq $ID1 $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "sbatch -o $LF/%A_%a.out -e $LF/%A_%a.err -J drtr-{}-$DATASET -a $ID1-$ID2 --nice=1000 bin/drtransformer.py data/$DATASET --id {}"

# if [ $ID1 -eq 0 ] && [ $ID2 -eq 0 ] # compute complete dataset at once
# then
#     sbatch -o $LF/%A.out -e $LF/%A.err -J drtr-$DATASET --nice=1000 bin/drtransformer.py data/$DATASET --all

# else # compute ID1 to ID2 (both inclusive), each as separate job
#     seq $ID1 $ID2 | xargs -I{} --max-procs=$NTHREADS bash -c "sbatch -o $LF/%A_%a.out -e $LF/%A_%a.err -J drtr-{}-$DATASET -a $ID1-$ID2 --nice=1000 bin/drtransformer.py data/$DATASET --id {}"
# fi
