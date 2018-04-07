% Produced by CVXGEN, 2018-04-03 18:09:48 -0400.
% CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2017 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: cvxsolve.m.
% Description: Solution file, via cvx, for use with sample.m.
function [vars, status] = cvxsolve(params, settings)
Hi = params.Hi;
cvx_begin
  % Caution: automatically generated by cvxgen. May be incorrect.
  variable Ui(15, 1);

  minimize(quad_form(Ui - Hi, eye(15)));
  subject to
    Ui >= 0;
    sum(Ui) == 1;
cvx_end
vars.Ui = Ui;
status.cvx_status = cvx_status;
% Provide a drop-in replacement for csolve.
status.optval = cvx_optval;
status.converged = strcmp(cvx_status, 'Solved');