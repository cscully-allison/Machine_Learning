#!/bin/bash

MAXROWS=$1
MAXGPUS=$2
MAXTRIALS=$3

echo "Numer of GPUs, Average Number of Iterations, Number Of Values Per GPU, Total Number of Values, Number Of Features, Number of Clusters, 8Avg. RSME, Total Time (ms), Time Per Iteration (ms)" > timingdata_multigpu_nrdc.csv

for((GPUS=1; GPUS <= $MAXGPUS; GPUS += 1))
do
    for((numVals=90; numVals < $MAXROWS; numVals *= 2))
    do
        for((trials = 0; trials < $MAXTRIALS; trials++))
        do
          srun --gres=gpu:$GPUS ./MAIN_PAR.py -i ../inputData/data.csv -o null -l 1 -k 15 -r .5  -t 50 -rw $numVals -c 11 -P .1 >> timingdata_multigpu_nrdc.csv
        done
    done
done
