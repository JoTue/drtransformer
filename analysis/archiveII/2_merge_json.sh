#!/bin/bash

DATASET=$1 #e.g. "tRNA"

# from: https://stackoverflow.com/questions/29636331/merging-json-files-using-a-bash-script
awk 'BEGIN{print "{"} FNR > 1 && last_file == FILENAME {print line} FNR == 1 {line = ""} FNR==1 && FNR != NR {printf ","} FNR > 1 {line = $0} {last_file = FILENAME} END{print "}"}' analysis_results/$DATASET/*.json > analysis_results/$DATASET/$DATASET.JSON && echo "json-files in" analysis_results/$DATASET "merged into" analysis_results/$DATASET/$DATASET.JSON
