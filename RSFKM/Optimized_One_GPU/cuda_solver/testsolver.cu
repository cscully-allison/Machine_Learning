/* Produced by CVXGEN, 2018-04-03 18:09:48 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.cu"

solver_scope do_solving(1);

void printVars(Vars vars);

#define NUMTESTS 0
int main(int argc, char **argv) {
  int num_iters;

  do_solving.set_defaults();
  do_solving.setup_indexing();
  do_solving.load_default_data();

  /* Solve problem instance for the record. */
  do_solving.settings.verbose = 1;
  num_iters = do_solving.solve();

  printVars(do_solving.vars);

  return 0;
}

void printVars(Vars vars){
	for(int x = 0; x < 15; x++){
		printf("Ui[%d]: %f \n", x, vars.Ui[x]);
	}
}
