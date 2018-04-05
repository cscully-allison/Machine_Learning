/* Produced by CVXGEN, 2018-04-03 18:09:48 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.h. */
/* Description: Header file with relevant definitions. */
#ifndef SOLVER_H
#define SOLVER_H
/* Uncomment the next line to remove all library dependencies. */
/*#define ZERO_LIBRARY_MODE */
#ifdef MATLAB_MEX_FILE
/* Matlab functions. MATLAB_MEX_FILE will be defined by the mex compiler. */
/* If you are not using the mex compiler, this functionality will not intrude, */
/* as it will be completely disabled at compile-time. */
#include "mex.h"
#else
#ifndef ZERO_LIBRARY_MODE
#include <stdio.h>
#endif
#endif
/* Space must be allocated somewhere (testsolver.c, csolve.c or your own */
/* program) for the global variables vars, params, work and settings. */
/* At the bottom of this file, they are externed. */
#ifndef ZERO_LIBRARY_MODE
#include <math.h>
#define pm(A, m, n) printmatrix(#A, A, m, n, 1)
#endif

typedef struct Params_t {
  double Hi[15];
} Params;
typedef struct Vars_t {
  double *Ui; /* 15 rows. */
} Vars;
typedef struct Workspace_t {
  double h[15];
  double s_inv[15];
  double s_inv_z[15];
  double b[1];
  double q[15];
  double rhs[46];
  double x[46];
  double *s;
  double *z;
  double *y;
  double lhs_aff[46];
  double lhs_cc[46];
  double buffer[46];
  double buffer2[46];
  double KKT[90];
  double L[45];
  double d[46];
  double v[46];
  double d_inv[46];
  double gap;
  double optval;
  double ineq_resid_squared;
  double eq_resid_squared;
  double block_33[1];
  /* Pre-op symbols. */
  double quad_640466485248[1];
  int converged;
} Workspace;
typedef struct Settings_t {
  double resid_tol;
  double eps;
  int max_iters;
  int refine_steps;
  int better_start;
  /* Better start obviates the need for s_init and z_init. */
  double s_init;
  double z_init;
  int verbose;
  /* Show extra details of the iterative refinement steps. */
  int verbose_refinement;
  int debug;
  /* For regularization. Minimum value of abs(D_ii) in the kkt D factor. */
  double kkt_reg;
} Settings;



//extern Vars vars;
//extern Params params;
//extern Workspace work;
//extern Settings settings;

/* Function definitions in ldl.c: */
__device__ void ldl_solve(double *target, double *var);
__device__ void ldl_factor(void);
__device__ double check_factorization(void);
__device__ void matrix_multiply(double *result, double *source);
__device__ double check_residual(double *target, double *multiplicand);
__device__ void fill_KKT(void);

/* Function definitions in matrix_support.c: */
__device__ void multbymA(double *lhs, double *rhs);
__device__ void multbymAT(double *lhs, double *rhs);
__device__ void multbymG(double *lhs, double *rhs);
__device__ void multbymGT(double *lhs, double *rhs);
__device__ void multbyP(double *lhs, double *rhs);
__device__ void fillq(Workspace *work, Params *params);
__device__ void fillh(Workspace *work);
__device__ void fillb(Workspace *work);
__device__ void pre_ops(Workspace *work, Params *params);

/* Function definitions in solver.c: */
__device__ double eval_gap(void);
__device__ void set_defaults(void);
__device__ void setup_pointers(void);
__device__ void setup_indexing(void);
__device__ void set_start(void);
__device__ double eval_objv(void);
__device__ void fillrhs_aff(void);
__device__ void fillrhs_cc(void);
__device__ void refine(double *target, double *var);
__device__ double calc_ineq_resid_squared(void);
__device__ double calc_eq_resid_squared(void);
__device__ void better_start(void);
__device__ void fillrhs_start(void);
__device__ long solve(void);

/* Function definitions in util.c: */
__device__ void tic(void);
__device__ float toc(void);
__device__ float tocq(void);
__device__ void printmatrix(char *name, double *A, int m, int n, int sparse);
__device__ double unif(double lower, double upper);
__device__ float ran1(long*idum, int reset);
__device__ float randn_internal(long *idum, int reset);
__device__ double randn(void);
__device__ void reset_rand(void);


#endif
