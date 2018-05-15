/* Produced by CVXGEN, 2018-04-03 18:09:48 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(1)-rhs[1]*(1)-rhs[2]*(1)-rhs[3]*(1)-rhs[4]*(1)-rhs[5]*(1)-rhs[6]*(1)-rhs[7]*(1)-rhs[8]*(1)-rhs[9]*(1)-rhs[10]*(1)-rhs[11]*(1)-rhs[12]*(1)-rhs[13]*(1)-rhs[14]*(1);
}
void multbymAT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(1);
  lhs[1] = -rhs[0]*(1);
  lhs[2] = -rhs[0]*(1);
  lhs[3] = -rhs[0]*(1);
  lhs[4] = -rhs[0]*(1);
  lhs[5] = -rhs[0]*(1);
  lhs[6] = -rhs[0]*(1);
  lhs[7] = -rhs[0]*(1);
  lhs[8] = -rhs[0]*(1);
  lhs[9] = -rhs[0]*(1);
  lhs[10] = -rhs[0]*(1);
  lhs[11] = -rhs[0]*(1);
  lhs[12] = -rhs[0]*(1);
  lhs[13] = -rhs[0]*(1);
  lhs[14] = -rhs[0]*(1);
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[1]*(-1);
  lhs[2] = -rhs[2]*(-1);
  lhs[3] = -rhs[3]*(-1);
  lhs[4] = -rhs[4]*(-1);
  lhs[5] = -rhs[5]*(-1);
  lhs[6] = -rhs[6]*(-1);
  lhs[7] = -rhs[7]*(-1);
  lhs[8] = -rhs[8]*(-1);
  lhs[9] = -rhs[9]*(-1);
  lhs[10] = -rhs[10]*(-1);
  lhs[11] = -rhs[11]*(-1);
  lhs[12] = -rhs[12]*(-1);
  lhs[13] = -rhs[13]*(-1);
  lhs[14] = -rhs[14]*(-1);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[1]*(-1);
  lhs[2] = -rhs[2]*(-1);
  lhs[3] = -rhs[3]*(-1);
  lhs[4] = -rhs[4]*(-1);
  lhs[5] = -rhs[5]*(-1);
  lhs[6] = -rhs[6]*(-1);
  lhs[7] = -rhs[7]*(-1);
  lhs[8] = -rhs[8]*(-1);
  lhs[9] = -rhs[9]*(-1);
  lhs[10] = -rhs[10]*(-1);
  lhs[11] = -rhs[11]*(-1);
  lhs[12] = -rhs[12]*(-1);
  lhs[13] = -rhs[13]*(-1);
  lhs[14] = -rhs[14]*(-1);
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2);
  lhs[1] = rhs[1]*(2);
  lhs[2] = rhs[2]*(2);
  lhs[3] = rhs[3]*(2);
  lhs[4] = rhs[4]*(2);
  lhs[5] = rhs[5]*(2);
  lhs[6] = rhs[6]*(2);
  lhs[7] = rhs[7]*(2);
  lhs[8] = rhs[8]*(2);
  lhs[9] = rhs[9]*(2);
  lhs[10] = rhs[10]*(2);
  lhs[11] = rhs[11]*(2);
  lhs[12] = rhs[12]*(2);
  lhs[13] = rhs[13]*(2);
  lhs[14] = rhs[14]*(2);
}
void fillq(void) {
  work.q[0] = -2*params.Hi[0];
  work.q[1] = -2*params.Hi[1];
  work.q[2] = -2*params.Hi[2];
  work.q[3] = -2*params.Hi[3];
  work.q[4] = -2*params.Hi[4];
  work.q[5] = -2*params.Hi[5];
  work.q[6] = -2*params.Hi[6];
  work.q[7] = -2*params.Hi[7];
  work.q[8] = -2*params.Hi[8];
  work.q[9] = -2*params.Hi[9];
  work.q[10] = -2*params.Hi[10];
  work.q[11] = -2*params.Hi[11];
  work.q[12] = -2*params.Hi[12];
  work.q[13] = -2*params.Hi[13];
  work.q[14] = -2*params.Hi[14];
}
void fillh(void) {
  work.h[0] = 0;
  work.h[1] = 0;
  work.h[2] = 0;
  work.h[3] = 0;
  work.h[4] = 0;
  work.h[5] = 0;
  work.h[6] = 0;
  work.h[7] = 0;
  work.h[8] = 0;
  work.h[9] = 0;
  work.h[10] = 0;
  work.h[11] = 0;
  work.h[12] = 0;
  work.h[13] = 0;
  work.h[14] = 0;
}
void fillb(void) {
  work.b[0] = 1;
}
void pre_ops(void) {
  work.quad_640466485248[0] = params.Hi[0]*params.Hi[0]+params.Hi[1]*params.Hi[1]+params.Hi[2]*params.Hi[2]+params.Hi[3]*params.Hi[3]+params.Hi[4]*params.Hi[4]+params.Hi[5]*params.Hi[5]+params.Hi[6]*params.Hi[6]+params.Hi[7]*params.Hi[7]+params.Hi[8]*params.Hi[8]+params.Hi[9]*params.Hi[9]+params.Hi[10]*params.Hi[10]+params.Hi[11]*params.Hi[11]+params.Hi[12]*params.Hi[12]+params.Hi[13]*params.Hi[13]+params.Hi[14]*params.Hi[14];
}
