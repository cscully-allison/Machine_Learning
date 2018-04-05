/* Produced by CVXGEN, 2018-04-03 18:09:48 -0400.  */
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
#define NUMTESTS 0
int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data();
  
  /* Solve problem instance for the record. */
  settings.verbose = 1;
  num_iters = solve();
  
  printVars(vars);

  
#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;
  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();
  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;
  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif
  return 0;
}

void printVars(Vars vars){
	for(int x = 0; x < 15; x++){
		printf("Ui[%d]: %f \n", x, vars.Ui[x]);
	}
}

void load_default_data(void) {
  params.Hi[0] = 0.20319161029830202;
  params.Hi[1] = 0.8325912904724193;
  params.Hi[2] = -0.8363810443482227;
  params.Hi[3] = 0.04331042079065206;
  params.Hi[4] = 1.5717878173906188;
  params.Hi[5] = 1.5851723557337523;
  params.Hi[6] = -1.497658758144655;
  params.Hi[7] = -1.171028487447253;
  params.Hi[8] = -1.7941311867966805;
  params.Hi[9] = -0.23676062539745413;
  params.Hi[10] = -1.8804951564857322;
  params.Hi[11] = -0.17266710242115568;
  params.Hi[12] = 0.596576190459043;
  params.Hi[13] = -0.8860508694080989;
  params.Hi[14] = 0.7050196079205251;
}
