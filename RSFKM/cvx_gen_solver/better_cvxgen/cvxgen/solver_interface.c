/* Produced by CVXGEN, 2018-04-11 10:54:25 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"

Vars vars;
Params params;
Workspace work;
Settings settings;

int main(int argc, char **argv) {
    int num_iters, NumRows, NumFeatures;
    double *H_Data;
    double RegParam;

    //Read in File

    //Tokenize input from File

    //store discrete Vars

    //allocate memory in H_Data

    //load H_Data with values


    set_defaults();
    setup_indexing();

    for(int row = 0, offset = 0; row < NumRows; row++, offset++){

        load_default_data(&H_Data[offset]);

        settings.verbose = 0;
        num_iters = solve();
        for(int i = 0; i < 15; i++){
            printf("Ui[%d] = %f \n", i, vars.Ui[i]);
        }
    }

    return 0;
}
void load_default_data(double *dataSubset) {
  params.Hi[0] = 14.704956;
  params.Hi[1] = 12.70112146;
  params.Hi[2] = 10.69518431;
  params.Hi[3] = 1.08427701;
  params.Hi[4] = 14.59945326;
  params.Hi[5] = 9.25445771;
  params.Hi[6] = 11.24174195;
  params.Hi[7] = 7.34465975;
  params.Hi[8] = 7.36713808;
  params.Hi[9] = 8.17758526;
  params.Hi[10] = 9.72512697;
  params.Hi[11] = 6.88102532;
  params.Hi[12] = 7.25736498;
  params.Hi[13] = 11.0101486;
  params.Hi[14] = 5.30128013;
}
