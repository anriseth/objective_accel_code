function [par, out] = ngmres_test_general(seednum)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('*********** NEW RUN ***********')
    %clear all;
    %------------------------------
    % problem description
    %------------------------------
    par.problem=3;
    % 3: standard quadratic function with diagonal matrix A and
    %    paraboloid coordinate transformation (problem B from the
    %    paper, arXiv:1106.4426v2)
    par.probPars{3}=[100 1e0];      % [size factor]

    %------------------------------
    % parameters
    %------------------------------
    % some general parameters
    % first select which methods to compare
    par.compareNGMRES_sdls_precond=1;       % run 1
    par.compareNGMRES_descent_precond=1;    % run 2
    par.compareNCG=1;                       % run 4
    par.compare_sdls_precond=0;             % run 5
    par.compareLBFGS=1;                     % run 7
    par.compareNGMRESO_sdls_precond=1;    % run 8
    par.compareNGMRESO_descent_precond=1; % run 9
                                          %
    par.figRunFirst=7; % which run to plot first in the figures (needs to be selected above)
    par.figStart=100;  % figure number
    par.initSeed=seednum;    % random seed
    s = RandStream('mt19937ar','Seed', par.initSeed);
    RandStream.setGlobalStream(s);
    par.precStep1=1e-4; % step in steepest descent preconditioner without line search
    par.precStep2=1e0;  % step factor in steepest descent preconditioner without line search
    par.maxfev=20;      % maximum function evaluations in all line searches
    par.restartfigure=false; % indicate restarts on the pictures
    par.epsi=1e-16;     % for plotting convergence to f*
                        % line styles
    par.fs{1}=':k';
    par.fs{2}='-k';
    par.fs{4}='-.k';
    par.fs{5}='--k';
    par.fs{7}='-vk';
    par.fs{8}=':ok';
    par.fs{9}='-ok';
    par.mispace{1} = 5;
    par.mispace{2} = 5;
    par.mispace{4} = 5;
    par.mispace{5} = 5;
    par.mispace{7} = 5;
    par.mispace{8} = 5;
    par.mispace{9} = 5;

    % ngmres parameters
    par.par_ngmres.w=20;       % maximum window size
    par.par_ngmres.maxIt=200;  % maximum number of iterations
    par.par_ngmres.relfTol=-1; % stopping tolerance on relative change in f
    par.par_ngmres.absgTol=-1; % stopping tolerance on the 2-norm of g
    par.par_ngmres.verbose=2;  % level of output to screen (0=no output, 1=summary output,
                               %       2=detailed output for each iteration)
    par.par_ngmres.logfev=true;% flag to include the history of the number of f and g
                               %       evaluations in the output (true or false)
    par.par_ngmres.logf=true;  % flag to include the history of functiona values in the output
                               %       (true or false)
    par.par_ngmres.logg=true;  % flag to include the history of the 2-norm of the gradient
    par.par_ngmres.logRestart=true; % flag to include the history of the restarts

    % ngmres line search parameters (for poblano_linesearch)
    par.par_ngmres.LineSearch_ftol=1.000000000000000e-04;
    par.par_ngmres.LineSearch_gtol=0.010000000000000;
    par.par_ngmres.LineSearch_initialstep=1;
    par.par_ngmres.LineSearch_maxfev=par.maxfev;
    par.par_ngmres.LineSearch_method='more-thuente';
    par.par_ngmres.LineSearch_stpmax=1.000000000000000e+15;
    par.par_ngmres.LineSearch_stpmin=1.000000000000000e-15;
    par.par_ngmres.LineSearch_xtol=1.000000000000000e-15;

    % ncg parameters (for comparing ncg with ngmres)
    par.par_ncg=ncg('defaults');
    par.par_ncg.MaxFuncEvals=50000;
    par.par_ncg.MaxIters=par.par_ngmres.maxIt;
    par.par_ncg.TraceFuncEvals=true;
    par.par_ncg.TraceFunc=true;
    par.par_ncg.TraceGradNorm=true;
    par.par_ncg.RelFuncTol=-1;
    par.par_ncg.StopTol=1e-300;
    par.par_ncg.Display='iter'; % iter, final or off
    par.par_ncg.RestartNW=false;
    par.par_ncg.Update='PR'; % 'PR', 'FR', 'HS' or 'SD'
    par.par_ncg.LineSearch_maxfev=par.maxfev; % default: 20
    par.par_ncg.LineSearch_ftol=1e-4; % default: 1e-4
    par.par_ncg.LineSearch_gtol=1e-2; % default: 1e-2
    par.par_ncg.LineSearch_initialstep=1;

    % lbfgs parameters (for comparing lbfgs with ngmres)
    par.par_lbfgs=lbfgs('defaults');
    par.par_lbfgs.M=5; % window size
    par.par_lbfgs.MaxFuncEvals=50000;
    par.par_lbfgs.MaxIters=par.par_ngmres.maxIt;
    par.par_lbfgs.TraceFuncEvals=true;
    par.par_lbfgs.TraceFunc=true;
    par.par_lbfgs.TraceGradNorm=true;
    par.par_lbfgs.RelFuncTol=-1;
    par.par_lbfgs.StopTol=1e-300;
    par.par_lbfgs.Display='iter'; % iter, final or off
    par.par_lbfgs.LineSearch_maxfev=par.maxfev; % default: 20
    par.par_lbfgs.LineSearch_ftol=1e-4; % default: 1e-4
    par.par_lbfgs.LineSearch_gtol=1e-2; % default: 1e-2
    par.par_lbfgs.LineSearch_initialstep=1;

    % steepest decent line search parameters (for using SD as a preconditioner)
    par.par_sdls_precond.LineSearch_ftol=1.000000000000000e-04;
    par.par_sdls_precond.LineSearch_gtol=0.010000000000000;
    par.par_sdls_precond.LineSearch_initialstep=1;
    par.par_sdls_precond.LineSearch_maxfev=par.maxfev;
    par.par_sdls_precond.LineSearch_method='more-thuente';
    par.par_sdls_precond.LineSearch_stpmax=1.000000000000000e+15;
    par.par_sdls_precond.LineSearch_stpmin=1.000000000000000e-15;
    par.par_sdls_precond.LineSearch_xtol=1.000000000000000e-15;

    %---------------------------------------------
    % set up the test problems
    %---------------------------------------------
    % the line search function for ngmres
    lineSearch_ngmres=@(fg,u0,f0,g0,d) poblano_linesearch(fg,u0,f0,g0,1.,d,par.par_ngmres);
    % the line search function for the sd preconditioner
    lineSearch_sdls_precond=@(fg,u0,f0,g0,d) poblano_linesearch(fg,u0,f0,g0,1.,d,par.par_sdls_precond);

    if par.problem==3
        % standard quadratic function with diagonal matrix A
        n=par.probPars{3}(1);
        kappa=par.probPars{3}(2);
        d=[1:n]';
        d(n)=d(n)*kappa;
        u0=rand(n,1); % generate the initial guess
        u_exact=ones(n,1);
        alpha=10.; % factor in the paraboloid coordinate transformation
        fg=@(u) func_problemB(u,d,u_exact,alpha); % set the objective function for N-GMRES
    end

    % the steepest descent preconditioner with line search
    M_sdls_precond=@(u,f,g) sd(u,f,g, ...
                               fg, ...                         % function that computes f and g
                               lineSearch_sdls_precond,...
                               par.par_sdls_precond);          % sd parameters
                                                               % the descent preconditioner with a fixed step
    M_descent_precond=@(u,f,g) descentDir(u,f,g, ...
                                          fg, ...                         % function that computes f and g
                                          par.precStep1,par.precStep2);

    %--------------------------------------
    % call the ngmres method
    %--------------------------------------
    if par.compareNGMRES_sdls_precond==1; % first N-GMRES, preconditioner: SD with line search
        disp('+++ start n-gmres with sdls preconditioner')
        out.out_ngmres_sdls_precond=ngmres(u0, ...            % initial guess
                                           fg, ...                % function that computes f and g
                                           M_sdls_precond, ...    % preconditioner function
                                           lineSearch_ngmres,...  % line search function
                                           par.par_ngmres);       % ngmres parameters
    end
    if par.compareNGMRES_descent_precond==1; % then N-GMRES, preconditioner: SD with small step
        disp('+++ start n-gmres with descent preconditioner')
        out.out_ngmres_descent_precond=ngmres(u0, ...         % initial guess
                                              fg, ...                % function that computes f and g
                                              M_descent_precond, ... % preconditioner function
                                              lineSearch_ngmres,...  % line search function
                                              par.par_ngmres);       % ngmres parameters
    end
    if par.compareNGMRESO_sdls_precond==1; % first N-GMRES, preconditioner: SD with line search
        disp('+++ start n-gmres-o with sdls preconditioner')
        out.out_ngmreso_sdls_precond=ngmres_o(u0, ...            % initial guess
                                              fg, ...                % function that computes f and g
                                              M_sdls_precond, ...    % preconditioner function
                                              lineSearch_ngmres,...  % line search function
                                              par.par_ngmres);       % ngmres parameters
    end
    if par.compareNGMRES_descent_precond==1; % then N-GMRES, preconditioner: SD with small step
        disp('+++ start n-gmreso with descent preconditioner')
        out.out_ngmreso_descent_precond=ngmres_o(u0, ...         % initial guess
                                                 fg, ...                % function that computes f and g
                                                 M_descent_precond, ... % preconditioner function
                                                 lineSearch_ngmres,...  % line search function
                                                 par.par_ngmres);       % ngmres parameters
    end


    %--------------------------------------
    % call ncg for comparison
    %--------------------------------------
    if par.compareNCG==1
        disp('+++ start n-cg')
        out.out_ncg=ncg(fg,u0,par.par_ncg);
        for i=2:size(out.out_ncg.TraceFuncEvals,2)
            out.out_ncg.TraceFuncEvals(i)=out.out_ncg.TraceFuncEvals(i)+out.out_ncg.TraceFuncEvals(i-1);
        end
    end

    %--------------------------------------
    % call lbfgs for comparison
    %--------------------------------------
    if par.compareLBFGS==1
        disp('+++ start lbfgs')
        out.out_lbfgs=lbfgs(fg,u0,par.par_lbfgs);
        for i=2:size(out.out_lbfgs.TraceFuncEvals,2)
            out.out_lbfgs.TraceFuncEvals(i)=out.out_lbfgs.TraceFuncEvals(i)+out.out_lbfgs.TraceFuncEvals(i-1);
        end
    end

    %--------------------------------------
    % compare with preconditioner only (steepest descent with line search)
    %--------------------------------------
    if par.compare_sdls_precond==1
        disp('+++ start sdls preconditioner only')
        u=u0;
        [f g]=fg(u);
        nfev=0;
        for k=1:par.par_ngmres.maxIt
            [u,f,g,fev]=M_sdls_precond(u,f,g);
            nfev=nfev+fev;
            out.out_sdls_precond.logf(k)=f;
            out.out_sdls_precond.logg(k)=norm(g);
            out.out_sdls_precond.logfev(k)=nfev;
        end
    end

    %-----------------------------------------
    % some figure output
    %--------------------------------------
    figure(par.figStart+2)
    util_ffigure(par,out)
    title('convergence towards the minimum value of f')

    figure(par.figStart+1)
    util_ffigure_evals(par,out)
    title('convergence towards the minimum value of f')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
