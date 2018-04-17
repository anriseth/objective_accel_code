function [par,out] = ngmres_test_tensor_CP(seednum, maxit, plotfigs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NGMRES_TEST_TENSOR_CP is a simple test script applying the N-GMRES
%optimization method with ALS preconditioning to a canonical tensor
%decomposition problem, and comparing with other optimization algorithms
%(NCG and ALS)
%
%   This function requires the Tensor toolbox (version 2.5) and the Poblano
%   toolbox. (http://www.sandia.gov/~tgkolda/TensorToolbox)
%
%   The nonlinear GMRES (N-GMRES) optimization method was proposed and analyzed
%   in the following references:
%      -Hans De Sterck, "A Nonlinear GMRES Optimization Algorithm for Canonical
%       Tensor Decomposition", SIAM J. Sci. Comput. 34, A1351-A1379, 2012.
%      -Hans De Sterck, "Steepest Descent Preconditioning for Nonlinear GMRES
%       Optimization", Numerical Linear Algebra with Applications 20, 453-471, 2013.
%
%   The preconditioner chosen for the N-GMRES algorithm is ALS here, and
%   poblano_linesearch from the Poblano toolbox is used for the line
%   searches.
%   NGMRES_TEST_TENSOR_CP calls the N-GMRES optimization method as
%   implemented in ngmres.m (via a call to CP_NGMRES).
%
%   The tensor test problem considered here is the problem from Tomasi and
%   Bro (2006) and Acar et al. (2011) for a ramdom dense tensor modified to
%   get specified collinearity between the columns.
%
%   See also NGMRES, CP_NGMRES, NGMRES_TEST_GENERAL, TENSOR, SPTENSOR,
%   KTENSOR, CP_OPT, NCG, POBLANO_LINESEARCH.
%
%by Hans De Sterck, September 2011
%edited by Asbjørn Nilsen Riseth, April 2018.
% - Modularised and added features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('*********** NEW RUN ***********')
    par = get_par(seednum,[],[],maxit);
    %% Compare different solvers than for the general problems
    par.compareNGMRESO_ALS = 1;
    par.compareNGMRES_ALS = 1;
    par.compareALS = 1;
    par.compareNGMRES_sdls = 0;
    par.compareNGMRESO_sdls = 0;
    par.compareNGMRES_sd = 0;
    par.compareNGMRESO_sd = 0;
    par.compareNCG = 1;
    par.compareLBFGS = 1;
    %% Get output results
    out = run_opts_tensor_CP(par);

    %-----------------------------------------
    % some figure output
    %--------------------------------------
    if plotfigs == true
        figure(par.figStart+2)
        util_ffigure(par,out)
        title('Tensor decomposition')

        figure(par.figStart+3)
        util_ffigure_evals(par,out)
        title('Tensor decomposition')
    end
end
