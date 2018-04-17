function par = get_par(seednum,probnum,n,maxit)
   %------------------------------
   % problem description
   %------------------------------
    par.problem=probnum;
    % 3: standard quadratic function with diagonal matrix A and
    %    paraboloid coordinate transformation (problem B from the
    %    paper, arXiv:1106.4426v2)
    par.probPars{1}=[n 0.0];     % problemA
    par.probPars{2}=[n 1.0 0.0]; % problemB
    par.probPars{3}=[n 0.0];     % problemC
    par.probPars{4}=[n 0.0];     % problemD
    par.probPars{5}=[n 0.0];     % problemE
    par.probPars{6}=[n 0.0];     % problemF
    par.probPars{7}=[n nan];     % problemG
    par.probPars{8}=[n 0.0];     % problemQuadratic
    par.probNames = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'Q'};

    %------------------------------
    % parameters
    %------------------------------
    % some general parameters
    % first select which methods to compare
    par.compareNGMRES_sdls=1;  % run 1
    par.compareNGMRES_sd=1;    % run 2
    par.compareNCG=1;          % run 4
    par.compare_sdls=0;        % run 5
    par.compareLBFGS=1;        % run 7
    par.compareNGMRESO_sdls=1; % run 8
    par.compareNGMRESO_sd=1;   % run 9

    % running parameters for the Tensor test problem
    par.compareNGMRESO_ALS=0;  % run 10
    par.compareNGMRES_ALS=0;   % run 11
    par.compareALS=0;          % run 12

    par.initSeed = seednum;

    par.figStart=100;        % figure number
    par.precStep1=1e-4;      % step in steepest descent preconditioner without line search
    par.precStep2=1e0;       % step factor in steepest descent preconditioner without line search
    par.maxfev=20;           % maximum function evaluations in all line searches
    par.restartfigure=false; % indicate restarts on the pictures
    par.epsi=1e-16;          % for plotting convergence to f*
                             % line styles
    par.fs{1}=':k';
    par.fs{2}='-k';
    par.fs{4}='-.k';
    par.fs{5}='--k';
    par.fs{7}='-vk';
    par.fs{8}=':ok';
    par.fs{9}='-ok';
    par.fs{10}=':k';  % overlaps with 1!
    par.fs{11}='--k'; % overlaps with 5!
    par.fs{12}=':ok'; % overlaps with 8!
    par.mispace{1} = 10;
    par.mispace{2} = 10;
    par.mispace{4} = 10;
    par.mispace{5} = 10;
    par.mispace{7} = 10;
    par.mispace{8} = 10;
    par.mispace{9} = 10;
    par.mispace{10} = 10;
    par.mispace{11} = 10;
    par.mispace{12} = 10;

    % ngmres parameters
    par.par_ngmres.epsi=1e-12;       % Regularization parameter for normal system
    par.par_ngmres.w=20;            % maximum window size
    par.par_ngmres.maxIt=maxit;      % maximum number of iterations
    par.par_ngmres.relfTol=-1;      % stopping tolerance on relative change in f
    par.par_ngmres.absgTol=1e-14;      % stopping tolerance on the
                                       % 2-norm of g (scaled by 1/n)
    par.par_ngmres.verbose=0;       % level of output to screen (0=no output, 1=summary output,
                                    %       2=detailed output for each iteration)
    par.par_ngmres.logfev=true;     % flag to include the history of the number of f and g
                                    %       evaluations in the output (true or false)
    par.par_ngmres.logf=true;       % flag to include the history of functiona values in the output
                                    %       (true or false)
    par.par_ngmres.logg=true;       % flag to include the history of the 2-norm of the gradient
    par.par_ngmres.logRestart=true; % flag to include the history of the restarts

    % ngmres line search parameters (for poblano_linesearch)
    par.par_ngmres.LineSearch_ftol=1.000000000000000e-04;
    par.par_ngmres.LineSearch_gtol=1e-1; % De Sterck paper default 1e-2
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
    par.par_ncg.RelFuncTol=par.par_ngmres.relfTol;
    par.par_ncg.StopTol=par.par_ngmres.absgTol;
    par.par_ncg.Display='off';                % iter, final or off
    par.par_ncg.RestartNW=false;
    par.par_ncg.Update='PR';                  % 'PR', 'FR', 'HS' or 'SD'
    par.par_ncg.LineSearch_maxfev=par.maxfev; % default: 20
    par.par_ncg.LineSearch_ftol=1e-4;         % default: 1e-4
    par.par_ncg.LineSearch_gtol=1e-1; % De Sterck paper default 1e-2
    par.par_ncg.LineSearch_initialstep=1;

    % lbfgs parameters (for comparing lbfgs with ngmres)
    par.par_lbfgs=lbfgs('defaults');
    par.par_lbfgs.M=5;                          % window size
    par.par_lbfgs.MaxFuncEvals=50000;
    par.par_lbfgs.MaxIters=par.par_ngmres.maxIt;
    par.par_lbfgs.TraceFuncEvals=true;
    par.par_lbfgs.TraceFunc=true;
    par.par_lbfgs.TraceGradNorm=true;
    par.par_lbfgs.RelFuncTol=par.par_ngmres.relfTol;
    par.par_lbfgs.StopTol=par.par_ngmres.absgTol;
    par.par_lbfgs.Display='off';                % iter, final or off
    par.par_lbfgs.LineSearch_maxfev=par.maxfev; % default: 20
    par.par_lbfgs.LineSearch_ftol=1e-4;         % default: 1e-4
    par.par_lbfgs.LineSearch_gtol=1e-1; % De Sterck paper default 1e-2
    par.par_lbfgs.LineSearch_initialstep=1;

    % steepest decent line search parameters (for using SD as a preconditioner)
    par.par_sdls.LineSearch_ftol=1.000000000000000e-04;
    par.par_sdls.LineSearch_gtol=1e-1; % De Sterck paper default 1e-2
    par.par_sdls.LineSearch_initialstep=1;
    par.par_sdls.LineSearch_maxfev=par.maxfev;
    par.par_sdls.LineSearch_method='more-thuente';
    par.par_sdls.LineSearch_stpmax=1.000000000000000e+15;
    par.par_sdls.LineSearch_stpmin=1.000000000000000e-15;
    par.par_sdls.LineSearch_xtol=1.000000000000000e-15;

    par.par_als.maxIt = par.par_ngmres.maxIt;
end
