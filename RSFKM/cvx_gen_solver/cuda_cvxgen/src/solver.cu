/* Produced by CVXGEN, 2018-04-03 18:09:48 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.c. */
/* Description: Main solver file. */
#include "solver.h"
__device__ double eval_gap(Workspace *work) {
  int i;
  double gap;
  gap = 0;
  for (i = 0; i < 15; i++)
    gap += work->z[i]*work->s[i];
  return gap;
}
__device__ void set_defaults(Settings *settings) {
  settings->resid_tol = 1e-6;
  settings->eps = 1e-4;
  settings->max_iters = 25;
  settings->refine_steps = 1;
  settings->s_init = 1;
  settings->z_init = 1;
  settings->debug = 0;
  settings->verbose = 1;
  settings->verbose_refinement = 0;
  settings->better_start = 1;
  settings->kkt_reg = 1e-7;
}
__device__ void setup_pointers(Workspace *work, Vars *vars) {
  work->y = work->x + 15;
  work->s = work->x + 16;
  work->z = work->x + 31;
  vars->Ui = work->x + 0;
}

// __device__ void setup_indexing(void) {
//   setup_pointers();
// }
__device__ void set_start(Workspace *work, Settings *settings) {
  int i;
  for (i = 0; i < 15; i++)
    work->x[i] = 0;
  for (i = 0; i < 1; i++)
    work->y[i] = 0;
  for (i = 0; i < 15; i++)
    work->s[i] = (work->h[i] > 0) ? work->h[i] : settings->s_init;
  for (i = 0; i < 15; i++)
    work->z[i] = settings->z_init;
}
__device__ double eval_objv(Workspace *work) {
  int i;
  double objv;
  /* Borrow space in work->rhs. */
  multbyP(work->rhs, work->x);
  objv = 0;
  for (i = 0; i < 15; i++)
    objv += work->x[i]*work->rhs[i];
  objv *= 0.5;
  for (i = 0; i < 15; i++)
    objv += work->q[i]*work->x[i];
  objv += work->quad_640466485248[0];
  return objv;
}
__device__ void fillrhs_aff(Workspace *work) {
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work->rhs;
  r2 = work->rhs + 15;
  r3 = work->rhs + 30;
  r4 = work->rhs + 45;
  /* r1 = -A^Ty - G^Tz - Px - q. */
  multbymAT(r1, work->y);
  multbymGT(work->buffer, work->z);
  for (i = 0; i < 15; i++)
    r1[i] += work->buffer[i];
  multbyP(work->buffer, work->x);
  for (i = 0; i < 15; i++)
    r1[i] -= work->buffer[i] + work->q[i];
  /* r2 = -z. */
  for (i = 0; i < 15; i++)
    r2[i] = -work->z[i];
  /* r3 = -Gx - s + h. */
  multbymG(r3, work->x);
  for (i = 0; i < 15; i++)
    r3[i] += -work->s[i] + work->h[i];
  /* r4 = -Ax + b. */
  multbymA(r4, work->x);
  for (i = 0; i < 1; i++)
    r4[i] += work->b[i];
}
__device__ void fillrhs_cc(Workspace *work) {
  int i;
  double *r2;
  double *ds_aff, *dz_aff;
  double mu;
  double alpha;
  double sigma;
  double smu;
  double minval;
  r2 = work->rhs + 15;
  ds_aff = work->lhs_aff + 15;
  dz_aff = work->lhs_aff + 30;
  mu = 0;
  for (i = 0; i < 15; i++)
    mu += work->s[i]*work->z[i];
  /* Don't finish calculating mu quite yet. */
  /* Find min(min(ds./s), min(dz./z)). */
  minval = 0;
  for (i = 0; i < 15; i++)
    if (ds_aff[i] < minval*work->s[i])
      minval = ds_aff[i]/work->s[i];
  for (i = 0; i < 15; i++)
    if (dz_aff[i] < minval*work->z[i])
      minval = dz_aff[i]/work->z[i];
  /* Find alpha. */
  if (-1 < minval)
      alpha = 1;
  else
      alpha = -1/minval;
  sigma = 0;
  for (i = 0; i < 15; i++)
    sigma += (work->s[i] + alpha*ds_aff[i])*
      (work->z[i] + alpha*dz_aff[i]);
  sigma /= mu;
  sigma = sigma*sigma*sigma;
  /* Finish calculating mu now. */
  mu *= 0.06666666666666667;
  smu = sigma*mu;
  /* Fill-in the rhs. */
  for (i = 0; i < 15; i++)
    work->rhs[i] = 0;
  for (i = 30; i < 46; i++)
    work->rhs[i] = 0;
  for (i = 0; i < 15; i++)
    r2[i] = work->s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
}
__device__ void refine(double *target, double *var, Workspace *work, Settings *settings) {
  int i, j;
  double *residual = work->buffer;
  double norm2;
  double *new_var = work->buffer2;
  for (j = 0; j < settings->refine_steps; j++) {
    norm2 = 0;
    matrix_multiply(residual, var, work);
    for (i = 0; i < 46; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
#ifndef ZERO_LIBRARY_MODE
    if (settings->verbose_refinement) {
      if (j == 0)
        printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
      else
        printf("After refinement we get squared norm %.6g.\n", norm2);
    }
#endif
    /* Solve to find new_var = KKT \ (target - A*var). */
    ldl_solve(residual, new_var, work, settings);
    /* Update var += new_var, or var += KKT \ (target - A*var). */
    for (i = 0; i < 46; i++) {
      var[i] -= new_var[i];
    }
  }
#ifndef ZERO_LIBRARY_MODE
  if (settings->verbose_refinement) {
    /* Check the residual once more, but only if we're reporting it, since */
    /* it's expensive. */
    norm2 = 0;
    matrix_multiply(residual, var, work);
    for (i = 0; i < 46; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
    if (j == 0)
      printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
    else
      printf("After refinement we get squared norm %.6g.\n", norm2);
  }
#endif
}
__device__ double calc_ineq_resid_squared(Workspace *work) {
  /* Calculates the norm ||-Gx - s + h||. */
  double norm2_squared;
  int i;
  /* Find -Gx. */
  multbymG(work->buffer, work->x);
  /* Add -s + h. */
  for (i = 0; i < 15; i++)
    work->buffer[i] += -work->s[i] + work->h[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 15; i++)
    norm2_squared += work->buffer[i]*work->buffer[i];
  return norm2_squared;
}
__device__ double calc_eq_resid_squared(Workspace *work) {
  /* Calculates the norm ||-Ax + b||. */
  double norm2_squared;
  int i;
  /* Find -Ax. */
  multbymA(work->buffer, work->x);
  /* Add +b. */
  for (i = 0; i < 1; i++)
    work->buffer[i] += work->b[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 1; i++)
    norm2_squared += work->buffer[i]*work->buffer[i];
  return norm2_squared;
}
__device__ void better_start(Workspace *work, Settings* settings) {
  /* Calculates a better starting point, using a similar approach to CVXOPT. */
  /* Not yet speed optimized. */
  int i;
  double *x, *s, *z, *y;
  double alpha;
  work->block_33[0] = -1;
  /* Make sure sinvz is 1 to make hijacked KKT system ok. */
  for (i = 0; i < 15; i++)
    work->s_inv_z[i] = 1;
  fill_KKT(work);
  ldl_factor(work, settings);
  fillrhs_start(work);
  /* Borrow work->lhs_aff for the solution. */
  ldl_solve(work->rhs, work->lhs_aff, work, settings);
  /* Don't do any refinement for now. Precision doesn't matter too much. */
  x = work->lhs_aff;
  s = work->lhs_aff + 15;
  z = work->lhs_aff + 30;
  y = work->lhs_aff + 45;
  /* Just set x and y as is. */
  for (i = 0; i < 15; i++)
    work->x[i] = x[i];
  for (i = 0; i < 1; i++)
    work->y[i] = y[i];
  /* Now complete the initialization. Start with s. */
  /* Must have alpha > max(z). */
  alpha = -1e99;
  for (i = 0; i < 15; i++)
    if (alpha < z[i])
      alpha = z[i];
  if (alpha < 0) {
    for (i = 0; i < 15; i++)
      work->s[i] = -z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 15; i++)
      work->s[i] = -z[i] + alpha;
  }
  /* Now initialize z. */
  /* Now must have alpha > max(-z). */
  alpha = -1e99;
  for (i = 0; i < 15; i++)
    if (alpha < -z[i])
      alpha = -z[i];
  if (alpha < 0) {
    for (i = 0; i < 15; i++)
      work->z[i] = z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 15; i++)
      work->z[i] = z[i] + alpha;
  }
}
__device__ void fillrhs_start(Workspace *work) {
  /* Fill rhs with (-q, 0, h, b). */
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work->rhs;
  r2 = work->rhs + 15;
  r3 = work->rhs + 30;
  r4 = work->rhs + 45;
  for (i = 0; i < 15; i++)
    r1[i] = -work->q[i];
  for (i = 0; i < 15; i++)
    r2[i] = 0;
  for (i = 0; i < 15; i++)
    r3[i] = work->h[i];
  for (i = 0; i < 1; i++)
    r4[i] = work->b[i];
}
__device__ long solve(Workspace *work, Params *params, Settings *settings, Vars* vars) {
  int i;
  int iter;
  double *dx, *ds, *dy, *dz;
  double minval;
  double alpha;
  work->converged = 0;
  setup_pointers(work, vars);
  pre_ops(work, params);
#ifndef ZERO_LIBRARY_MODE
  if (settings->verbose)
    printf("iter     objv        gap       |Ax-b|    |Gx+s-h|    step\n");
#endif
  fillq(work, params);
  fillh(work);
  fillb(work);
  if (settings->better_start)
    better_start(work, settings);
  else
    set_start(work, settings);
  for (iter = 0; iter < settings->max_iters; iter++) {
    for (i = 0; i < 15; i++) {
      work->s_inv[i] = 1.0 / work->s[i];
      work->s_inv_z[i] = work->s_inv[i]*work->z[i];
    }
    work->block_33[0] = 0;
    fill_KKT(work);
    ldl_factor(work, settings);
    /* Affine scaling directions. */
    fillrhs_aff(work);
    ldl_solve(work->rhs, work->lhs_aff, work, settings);
    refine(work->rhs, work->lhs_aff, work, settings);
    /* Centering plus corrector directions. */
    fillrhs_cc(work);
    ldl_solve(work->rhs, work->lhs_aff, work, settings);
    refine(work->rhs, work->lhs_aff, work, settings);
    /* Add the two together and store in aff. */
    for (i = 0; i < 46; i++)
      work->lhs_aff[i] += work->lhs_cc[i];
    /* Rename aff to reflect its new meaning. */
    dx = work->lhs_aff;
    ds = work->lhs_aff + 15;
    dz = work->lhs_aff + 30;
    dy = work->lhs_aff + 45;
    /* Find min(min(ds./s), min(dz./z)). */
    minval = 0;
    for (i = 0; i < 15; i++)
      if (ds[i] < minval*work->s[i])
        minval = ds[i]/work->s[i];
    for (i = 0; i < 15; i++)
      if (dz[i] < minval*work->z[i])
        minval = dz[i]/work->z[i];
    /* Find alpha. */
    if (-0.99 < minval)
      alpha = 1;
    else
      alpha = -0.99/minval;
    /* Update the primal and dual variables. */
    for (i = 0; i < 15; i++)
      work->x[i] += alpha*dx[i];
    for (i = 0; i < 15; i++)
      work->s[i] += alpha*ds[i];
    for (i = 0; i < 15; i++)
      work->z[i] += alpha*dz[i];
    for (i = 0; i < 1; i++)
      work->y[i] += alpha*dy[i];
    work->gap = eval_gap(work);
    work->eq_resid_squared = calc_eq_resid_squared(work);
    work->ineq_resid_squared = calc_ineq_resid_squared(work);
#ifndef ZERO_LIBRARY_MODE
    if (settings->verbose) {
      work->optval = eval_objv(work);
      printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter+1, work->optval, work->gap, sqrt(work->eq_resid_squared),
          sqrt(work->ineq_resid_squared), alpha);
    }
#endif
    /* Test termination conditions. Requires optimality, and satisfied */
    /* constraints. */
    if (   (work->gap < settings->eps)
        && (work->eq_resid_squared <= settings->resid_tol*settings->resid_tol)
        && (work->ineq_resid_squared <= settings->resid_tol*settings->resid_tol)
       ) {
      work->converged = 1;
      work->optval = eval_objv(work);
      return iter+1;
    }
  }
  return iter;
}
