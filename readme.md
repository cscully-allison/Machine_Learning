#About
Author: Connor Scully-Allison
Date: 05/15/2018
This is a repository for all code and data related to the multi-GPU implementation of a robust and sparse fuzzy-k means algorithim.

#Main Directories
1. inputData - This directory holds all data which was used for testing at various points of development. Timings were collected using data.csv.
2. Seqential - This directory contains the executable sequential code of RSFKM used to collect baseline timing data. The subdirectory cvxgen containes generated optimization code.
3. Multiple-GPU - This diectory contains the executable GPU optimized code for iRSFKM. The subdirectory cvx gen contains CPU adapatations of generated opimization code and all the source code for all cuda kernels called from RSFKM.py. 
4. legacy - Contains various iterations of RSFKM which were not used for data colletion.


#Running the program

##Setup
From RSFKM/ run command

```
export CUDA_HOME=/usr/local/cuda-9.0
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64 
PATH=${CUDA_HOME}/bin:${PATH}
export PATH

source PyEnv/bin/activate
```

##Seqential

From the Sequential directory, optimized sequential code can be run with the following command:

```
./MAIN.py -i ../inputData/data.csv -o null -l 1 -k 15 -r .5 -t 100 -rw 1000 -c 11 -g 1
```

* -i: input data
* -o: output of graph to show clustering, does not work on cubix box
* -l: indicates if there is a leading column of time stamps which needs to be stripped away for clustering
* -k: Number of clusters that our data will be clustered into; with the optimized CPU code, this is must be set to 15
* -r: Regukating parameter used to enforce sparseness; can be tuned to produce different clustering results
* -t: Threshold parameter used to reduce influence of outlers on centroid updates; can be tuned to produce different clustering resultsl sufficently low values <50 will result in errors however
* -rw: the number of rows we are clustering
* -c: the number of features we are using in our clustering; with data.csv this value can be up to 30


##Multi-Gpu

From the Multiple_GPU directory, GPU code can be run with the following command:

```
srun --gres=gpu:1 ./MAIN_PAR.py -i ../inputData/data.csv -o null -l 1 -k 15 -r .5  -t 50 -rw 1000 -c 11 -P .1
```

* --gres=gpu:1 : Sets the number of GPUs to be used. On cubix this can be any number between 1 and 8.
* -i: input data
* -o: output of graph to show clustering, does not work on cubix box
* -l: indicates if there is a leading column of time stamps which needs to be stripped away for clustering
* -k: Number of clusters that our data will be clustered into; with the optimized CPU code, this is must be set to 15
* -r: Regukating parameter used to enforce sparseness; can be tuned to produce different clustering results
* -t: Threshold parameter used to reduce influence of outlers on centroid updates; can be tuned to produce different clustering resultsl sufficently low values < 50 will result in errors however
* -rw: the number of rows we are clustering
* -c: the number of features we are using in our clustering; with data.csv this value can be up to 30
* -P: this argument denotes the percentage of values we are removing to test the accuracy of our imputation.


# Some Resources

http://stanford.edu/~cpiech/cs221/handouts/kmeans.html

https://en.wikipedia.org/wiki/Matrix_norm#Frobenius_norm
