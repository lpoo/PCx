%
% MATLAB-PCx interface v0.1 BETA
%
% Michael Wagner
% Department of Mathematics and Statistics
% Old Dominion University
% August 2000
%
% revised at CCHMC, Cincinnati, May 2003
%
%
%
%
% Syntax:
%
% [x,f,exitflag,output,lambda] = PCx(c, A_ineq, b_ineq, A_eq, b_eq, lb, ub, x0, opts)
%
% All but the last argument are mandatory. Use empty arrays if
% appropriate. The options-structure replaces the .specs-file.
%
% Input: c        - n-dimensional real array
%        A_ineq   - m2-by-n-dimensional SPARSE real array
%        b_ineq   - m2-dimensional real array
%        A_eq     - m1-by-n-dimensional SPARSE real array
%        b_eq     - m1-dimensional real array
%        lb       - n-dimensional real array, possibly -inf
%        ub       - n-dimensional real array, possibly +inf
%        x0       - starting point (IGNORED, for compatibility purposes only)
%        opts     - structure that can contain elements with names out of
%                  the following list: "max", "presolve", "opttol", 
%	               "prifeastol", "dualfeastol", "cachesize", 
%                  "unrollinglevel", "iterationlimit", "centerexp", 
%                  "refinement", "stepfactor", "scaling",
%                  "hocorrections", "maxcorrections". The corresponding
%                  values sould make sense, please see the user guide for
%                  a description. (Note: for boolean options we use the
%                  C-convention for  0 = true, 1 (or any other nonzero) =
%                  false). 
%
%                  The options "Display", "Diagnostics", "TolFun", and "MaxIter"
%                  are also recognized, so that overall the syntax is
%                  (almost) entirely compatible with that of 'linprog'
%                  from the MATLAB optimization toolbox.
%
% Output: primal-dual solution pair for the LP problem
%
%                         min c'*x
%                 s.t.
%                         A_eq * x  = b_eq
%                       A_ineq * x <= b_ineq
%                          lb <= x <= ub
%
%  
%  where x        - primal solution
%        f        - optimal value of LP 
%        exitflag - 1 for successful run, 0 or -1 for error
%        output   - structure with the number of iterations taken in
%                   output.iterations (and other, redundant field members).
%        lambda   - structure with the set of Lagrangian multipliers
%                   at the solution: 
%                        lambda.ineqlin for the linear inequalities A_ineq, 
%                        lambda.eqlin for the linear equalities A_eq,
%                        lambda.lower for LB, 
%                    and lambda.upper for UB
%  
% The user should make use of 'inf' for infinite bounds, avoid using just
% very large numbers. This can make the difference between a failure and
% success. 
%
