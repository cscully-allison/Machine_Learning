#!/bin/bash

MAXROWS=$1
MAXTRIALS=$2

echo "Number of Iterations, Number Of Values, Number Of Dimensions, Number of Clusters, Total Time (ms), Time Per Iteration (ms)" > timingdata_variablenumvals_nrdc.csv

for((numVals=64; numVals < $MAXROWS; numVals *= 2))
do
    for((trials = 0; trials < $MAXTRIALS; trials++))
    do
        ./MAIN.py -i ../inputData/data.csv -o null -l 1 -k 15 -r 1 -t 100 -rw $numVals -c 2 -g 0 >> timingdata_variablenumvals_nrdc.csv
    done

    if (($numVals == 100)); then
        numVals=0
    fi
done

#
# echo "Number of Iterations, Number Of Values, Number Of Dimensions, Number of Clusters, Total Time (ms), Time Per Iteration (ms)" > timingdata_variablenumdimensions.csv
#
# for((numDims=2; numDims < $MAXCOLS; numDims += 2))
# do
#
# done
