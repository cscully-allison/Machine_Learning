/* Produced by CVXGEN, 2018-04-11 10:54:25 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Vars vars;
Params params;
Workspace work;
Settings settings;

int main(int argc, char **argv) {
    int num_iters, NumRows, NumFeatures, filesize;
    char *buffer;
    double *H_Data;
    double *MM_Data;
    double RegParam;
    FILE *fp;

    //Read in File
    fp = fopen("interface_file.csv", "r");
    fseek(fp, 0L, SEEK_END);
    filesize = ftell(fp);
    rewind(fp);

    buffer = malloc(filesize*sizeof(char));
    fgets(buffer, filesize, fp);

    // //store discrete Vars
    NumRows = atoi( strtok(buffer, ",") );
    NumFeatures = atoi( strtok(NULL, ",") );
    RegParam = (double) atof( strtok(NULL, ",") );


    //allocate memory in H_Data
    H_Data = malloc( NumRows * NumFeatures * sizeof(double) );
    MM_Data = malloc( NumRows * NumFeatures * sizeof(double) );


    //load H_Data with values
    for(int i = 0; i < NumRows*NumFeatures; i++){
        H_Data[i] = (double) atof( strtok(NULL, ",") );
    }




    set_defaults();
    setup_indexing();

    for(int row = 0, offset = 0; row < NumRows; row++, offset += 15){

        load_default_data(&H_Data[offset], NumFeatures, RegParam);

        settings.verbose = 0;
        num_iters = solve();

        for(int i = 0; i < NumFeatures; i++){
            MM_Data[offset + i] = vars.Ui[i];

            printf("%f,", MM_Data[offset+i]);
        }




    }

    return 0;
}

void load_default_data(double *dataSubset, int NumFeatures, double RegParam) {

    double multiplicand = (-1/(2*RegParam));
    for(int i = 0; i < NumFeatures; i++){
        params.Hi[i] = dataSubset[i] * multiplicand;
    }

}
