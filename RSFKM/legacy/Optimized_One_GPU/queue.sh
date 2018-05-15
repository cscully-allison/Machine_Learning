#!/bin/bash

MAXROWS=$1
MAXGPUS=$2
MAXTRIALS=$3

echo "Number of Iterations, Number Of Values, Number Of Dimensions, Number of Clusters, Total Time (ms), Time Per Iteration (ms), Num Threads (from user), Total Num Threads (per Block) " > timingdata_multigpu_nrdc.csv

for((GPUS=1; GPUS < $MAXGPUS; GPUS += 1))
do
    for((numVals=64; numVals < $MAXROWS; numVals *= 2))
    do
        for((trials = 0; trials < $MAXTRIALS; trials++))
        do
          srun --gres=gpu:$GPUS ./MAIN_PAR.py -i ../inputData/data.csv -o null -l 1 -k 15 -r .5  -t 50 -rw 6000 -c 11 -P .01 >> timingdata_multigpu_nrdc.csv
        done

        if (($numVals == 100)); then
            numVals=0
        fi
    done
done
