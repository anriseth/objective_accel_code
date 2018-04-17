function [par,out] = ngmres_test_general(seednum,probnum,n,maxit,plotfigs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NGMRES_TEST_GENERAL is a simple test script applying the N-GMRES
%optimization method to a simple nonlinear optimization test problem, and
%comparing with other optimization algorithms (NCG, LBFGS, steepest descent)
%
%   This function requires the Poblano toolbox (for the linesearch and for
%   comparison with NCG and LBFGS). (software.sandia.gov/trac/poblano)
%
%   The nonlinear GMRES (N-GMRES) optimization method was proposed and analyzed
%   in the following references:
%      -Hans De Sterck, "A Nonlinear GMRES Optimization Algorithm for Canonical
%       Tensor Decomposition", SIAM J. Sci. Comput. 34, A1351-A1379, 2012.
%      -Hans De Sterck, "Steepest Descent Preconditioning for Nonlinear GMRES
%       Optimization", Numerical Linear Algebra with Applications 20, 453-471, 2013.
%
%   NGMRES_TEST_GENERAL calls the N-GMRES optimization method as
%   implemented in ngmres.m.
%
%   The test problem considered here is problem B from "Steepest Descent
%   Preconditioning for Nonlinear GMRES Optimization", Numerical Linear Algebra
%   with Applications 20, 453-471, 2013. As in that paper, we compare two
%   versions of the steepest descent preconditioner in this test script, and
%   poblano_linesearch from the Poblano toolbox is used for the line
%   searches.
%
%   NGMRES_TEST_GENERAL can easily be adapted to test N-GMRES for your own
%   unconstrained nonlinear optimization problem: just replace the function
%   'func_problemB' below with any other function that computes the
%   function value and gradient of an objective function.
%
%   See also NGMRES, POBLANO_LINESEARCH.
%
%by Hans De Sterck, September 2011
%
% Extensions by Asbjørn Nilsen Riseth, August 2017
% - Modularised and added features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('*********** NEW RUN ***********')
    par = get_par(seednum,probnum,n,maxit);
    out = run_opts(par);

    %-----------------------------------------
    % some figure output
    %--------------------------------------
    if plotfigs == true
        figure(par.figStart)
        util_ffigure(par,out)
        title(sprintf('Problem %s, seed number %d', par.probNames{probnum}, ...
                      seednum))

        figure(par.figStart+1)
        util_ffigure_evals(par,out)
        title(sprintf('Problem %s, seed number %d', par.probNames{probnum}, ...
                      seednum))
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
