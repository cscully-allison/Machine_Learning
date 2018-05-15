#!/bin/bash

MAXROWS=$1
MAXTRIALS=$2

echo "Number of Iterations, Number Of Values, Number Of Dimensions, Number of Clusters, Total Time (ms), Time Per Iteration (ms)" > timingdata_variablenumvals_nrdc.csv

for((numVals=90; numVals < $MAXROWS; numVals *= 2))
do
    for((trials = 0; trials < $MAXTRIALS; trials++))
    do
        ./MAIN.py -i ../inputData/data.csv -o null -l 1 -k 15 -r .5 -t 100 -rw $numVals -c 11 -g 1 >> timingdata_variablenumvals_nrdc.csv
    done
done

