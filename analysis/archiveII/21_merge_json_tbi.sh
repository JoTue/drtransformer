#!/bin/bash

DATASET=$1 #e.g. "tRNA"
DIRECTORY=$2

# from: https://stackoverflow.com/questions/29636331/merging-json-files-using-a-bash-script
awk 'BEGIN{print "{"} FNR > 1 && last_file == FILENAME {print line} FNR == 1 {line = ""} FNR==1 && FNR != NR {printf ","} FNR > 1 {line = $0} {last_file = FILENAME} END{print "}"}' ../$DIRECTORY/analysis_results/$DATASET/*.json > ../$DIRECTORY/analysis_results/$DATASET/$DATASET.JSON && echo "json-files in" ../$DIRECTORY/analysis_results/$DATASET "merged into" ../$DIRECTORY/analysis_results/$DATASET/$DATASET.JSON
