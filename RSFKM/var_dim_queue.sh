#!/bin/bash

MAXCOLS=$1
COLTAKE=$2
MAXTRIALS=$3

echo "Number of Iterations, Number Of Values, Number Of Dimensions, Number of Clusters, Total Time (ms), Time Per Iteration (ms)" > timingdata_variabledims.csv

for((numDims=2; numDims < $MAXCOLS; numDims += $COLTAKE))
do
    for((trials = 0; trials < $MAXTRIALS; trials++))
    do
        ./MAIN.py -i inputData/data.csv -o null -l 1 -k 4 -r 0.5 -t 100 -rw 1000 -c $numDims >> timingdata_variabledims.csv
    done

    if (($numDims == 1)); then
        numDims=0
    fi
done

#
# echo "Number of Iterations, Number Of Values, Number Of Dimensions, Number of Clusters, Total Time (ms), Time Per Iteration (ms)" > timingdata_variablenumdimensions.csv
#
# for((numDims=2; numDims < $MAXCOLS; numDims += 2))
# do
#
# done
