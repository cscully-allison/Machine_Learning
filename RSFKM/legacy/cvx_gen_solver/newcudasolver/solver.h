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
#define ZERO_LIBRARY_MODE

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

struct solver_scope{
    Vars vars;
    Params params;
    Workspace work;
    Settings settings;
    int id;

    solver_scope(int _id);

    /* Function definitions in ldl.c: */
    void ldl_solve(double *target, double *var);
    void ldl_factor(void);
    double check_factorization(void);
    void matrix_multiply(double *result, double *source);
    double check_residual(double *target, double *multiplicand);
    void fill_KKT(void);

    /* Function definitions in matrix_support.c: */
    void multbymA(double *lhs, double *rhs);
    void multbymAT(double *lhs, double *rhs);
    void multbymG(double *lhs, double *rhs);
    void multbymGT(double *lhs, double *rhs);
    void multbyP(double *lhs, double *rhs);
    void fillq(void);
    void fillh(void);
    void fillb(void);
    void pre_ops(void);

    /* Function definitions in solver.c: */
    double eval_gap(void);
    void set_defaults(void);
    void setup_pointers(void);
    void setup_indexing(void);
    void set_start(void);
    double eval_objv(void);
    void fillrhs_aff(void);
    void fillrhs_cc(void);
    void refine(double *target, double *var);
    double calc_ineq_resid_squared(void);
    double calc_eq_resid_squared(void);
    void better_start(void);
    void fillrhs_start(void);
    long solve(void);

    /* Function definitions in testsolver.c: */
    int main(int argc, char **argv);
    void load_default_data(void);
}
