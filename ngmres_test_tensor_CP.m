function ngmres_test_tensor_CP
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('*********** NEW RUN ***********')
    clear all;
    %------------------------------
    % parameters
    %------------------------------
    % some general parameters
    compareNGMRESO=1;
    compareNCG=1;
    compareALS=1;
    figStart=0;
    initSeed=0;
    s = RandStream('mt19937ar','Seed', initSeed);
    RandStream.setGlobalStream(s);

    % ngmres parameters
    par.w=20;      % maximum window size
    par.maxIt=200; % maximum number of iterations
    par.relfTol=-1;% stopping tolerance on relative change in f
    par.absgTol=-1;% stopping tolerance on the 2-norm of g
    par.verbose=2; % level of output to screen (0=no output, 1=summary output,
                   %       2=detailed output for each iteration)
    par.logfev=true; % flag to include the history of the number of f and g
                     %       evaluations in the output (true or false)
    par.logf=true; % flag to include the history of functiona values in the output
                   %       (true or false)
    par.logg=true; % flag to include the history of the 2-norm of the gradient
    par.logRestart=true; % flag to include the history of the restarts

    % ngmres line search parameters (for poblano_linesearch)
    par_ls.LineSearch_ftol=1.000000000000000e-04;
    par_ls.LineSearch_gtol=0.010000000000000;
    par_ls.LineSearch_initialstep=1;
    par_ls.LineSearch_maxfev=20;
    par_ls.LineSearch_method='more-thuente';
    par_ls.LineSearch_stpmax=1.000000000000000e+15;
    par_ls.LineSearch_stpmin=1.000000000000000e-15;
    par_ls.LineSearch_xtol=1.000000000000000e-15;

    % ncg parameters (for comparing ncg with ngmres)
    options=ncg('defaults');
    options.MaxFuncEvals=2000;
    options.MaxIters=par.maxIt;
    options.TraceFunc=true;
    options.TraceGradNorm=true;
    options.RelFuncTol=-1;
    options.StopTol=1e-99;
    options.Display='iter'; % iter, final or off
    options.RestartNW=false;
    options.Update='PR'; % 'PR', 'FR', 'HS' or 'SD'
    options.LineSearch_maxfev=20; % default: 20
    options.LineSearch_ftol=1e-4; % default: 1e-4
    options.LineSearch_gtol=1e-2; % default: 1e-2
    options.LineSearch_initialstep=1;

    %------------------------------
    % create a test tensor
    % (this is the test from Tomasi and Bro (2006) and Acar et al. (2011) for
    % a ramdom dense tensor modified to get specified collinearity between the
    % columns)
    %------------------------------
    s=50; % size of each mode
    r=3; % number of rank-one components in the original noise-free random tensor
    l1=1; % magnitude of the first type of noise
    l2=1; % magnitude of the second type of noise
    c=0.9; % collinearity
    pars=[s c r l1 l2];
    T=can_createTensorDense(pars);

    %--------------------------------------
    % do some pre-processing
    %--------------------------------------
    nT2=norm(T)^2;

    %--------------------------------------
    % generate the random initial guess
    %--------------------------------------
    u0=rand(r*sum(size(T)),1);

    %--------------------------------------
    % call cp_ngmres
    %--------------------------------------
    out=cp_ngmres(T,r,tt_cp_vec_to_fac(u0,T),par,par_ls);

    if compareNGMRESO==1
        out_o=cp_ngmres_o(T,r,tt_cp_vec_to_fac(u0,T),par,par_ls);
    end

    %--------------------------------------
    % call ncg for comparison
    %--------------------------------------
    if compareNCG==1
        disp('+++ start n-cg directly')
        out_ncg=ncg(@(u) tt_cp_fun(u,T,nT2),u0,options);
    end

    %--------------------------------------
    % compare with ALS
    %--------------------------------------
    if compareALS==1
        disp('+++ start ALS')
        u=u0;
        [f g]=tt_cp_fun(u,T,nT2);
        for k=1:par.maxIt
            [u,f,g,fev]=can_ALSu(T,r,u,f,g,@(u) tt_cp_fun(u,T,nT2));
            fALS(k)=f;
            gALS(k)=norm(g);
        end
    end

    %-----------------------------------------
    % some figure output
    %--------------------------------------
    if par.logg
        figure(figStart+1)
        semilogy(out.logg,'-+')
        title('convergence of the norm of the gradient')
        if par.logRestart
            loggMod=out.logg;
            loggMod(out.logRestart==0)=NaN;
            hold on
            semilogy(loggMod,'+r')
            hold off
        end
        if compareNGMRESO==1
            hold on
            semilogy(out_o.logg,'-o')
            hold off
        end
        if par.logRestart
            loggMod=out_o.logg;
            loggMod(out_o.logRestart==0)=NaN;
            hold on
            semilogy(loggMod,'+r')
            hold off
        end
        if compareNCG==1
            hold on
            semilogy(out_ncg.TraceGradNorm(2:end),'-*')
            hold off
        end
        if compareALS==1
            hold on
            semilogy(gALS,'-')
            hold off
        end
        if par.logRestart
            legend('n-gmres','restarts','n-gmres-o','restarts-o','n-cg','ALS')
        else
            legend('n-gmres','n-gmres-o','n-cg','ALS')
        end
        xlabel('iterations')
    end

    if par.logf
        figure(figStart+2)

        minval = min(min(out.logf),min(out_o.logf));
        minval = min(minval,min(out_ncg.TraceFunc(2:end)));
        minval = min(minval,min(fALS));
        semilogy(out.logf-minval,'-+')
        title('convergence towards the minimum value of f')
        if par.logRestart
            logfMod=out.logf-minval;
            logfMod(out.logRestart==0)=NaN;
            hold on
            semilogy(logfMod,'+r')
            hold off
        end
        if compareNGMRESO==1
            hold on
            semilogy(out_o.logf-minval,'-o')
            hold off
        end
        if par.logRestart
            logfMod=out_o.logf-minval;
            logfMod(out_o.logRestart==0)=NaN;
            hold on
            semilogy(logfMod,'+r')
            hold off
        end

        if compareNCG==1
            hold on
            semilogy(out_ncg.TraceFunc(2:end)-minval,'-*')
            hold off
        end
        if compareALS==1
            hold on
            semilogy(fALS-minval,'-')
            hold off
        end
        if par.logRestart
            legend('n-gmres','restarts','n-gmres-o','restarts-o','n-cg','ALS')
        else
            legend('n-gmres','n-gmres-o','n-cg','ALS')
        end
        xlabel('iterations')
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,f,g,fev]=can_ALSu(T,r,u,f,g,fg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function performs one ALS iteration (it uses a modification of the
%ALS routine in the tensor toolbox (only one iteration, and normalization
%is different))
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract number of dimensions and norm of T.
    nModes = ndims(T);
    dims=size(T);

    % convert u to a ktensor
    A = tt_cp_vec_to_fac(u,T);
    A=ktensor(A);

    A=can_ALS(T,r,A);

    % convert A to a vector
    u = tt_fac_to_vec(A.U);

    % assume ALS is equivalent to 2 f and g evaluations (in reality it is a
    % bit less expensive I think)
    fev=2;

    [f g]=fg(u);
    fev=fev+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A=can_ALS(T,r,A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function performs one ALS iteration
%this is a modification of the ALS routine in the tensor toolbox (only one
%iteration, and normalization is different)
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract number of dimensions and norm of T.
    nModes = ndims(T);

    % if there is a lambda factor in A, put it into the first mode
    lambda=A.lambda;
    A.U{1} = A.U{1} * spdiags(lambda,0,r,r);
    A.lambda=ones(r,1);

    % Iterate over all modes of the tensor
    for n = 1:nModes
        % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
        Unew = mttkrp(T,A.U,n);

        % Compute the matrix of coefficients for linear system
        Y = ones(r,r);
        for i = [1:n-1,n+1:nModes]
            Y = Y .* (A.U{i}'*A.U{i});
        end

        Unew = (Y \ Unew')'; %<- Line from TTB 2.2.

        if issparse(Unew)
            A.U{n} = full(Unew);   % for the case R=1
        else
            A.U{n} = Unew;
        end
    end

    A = ktensor(A.U); % default: lambda=1

    A=arrange(A); % normalize and introduce lambda, and make sure the components appear in order of lambda

    % now distribute lambda evenly over all modes
    facMat=spdiags(nthroot(A.lambda,nModes),0,r,r);
    for n=1:nModes
        A.U{n} = A.U{n} * facMat;
    end

    A = ktensor(A.U); % lambda=1 by default
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Zdprime=can_createTensorDense(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function creates a test tensor for CP decomposition
%(this is the test from Tomasi and Bro (2006) and Acar et al. (2011) for
%a ramdom dense tensor modified to get specified collinearity between the
%columns)
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sizet=params(1);
    c=params(2);
    rTrue=params(3);
    l1=params(4);
    l2=params(5);

    % first generate K
    K=ones(rTrue,rTrue)*c;
    for k=1:rTrue
        K(k,k)=1;
    end
    K;

    % get C as the cholesky factor of K
    C=chol(K);

    for n=1:3
        % generate a random matrix
        M=randn(sizet,rTrue);
        % ortho-normalize the columns of M, gives Q
        [Q,R] = qr(M,0);
        U{n}=Q*C;
    end

    Z=full(ktensor(U));

    % generate two random tensors
    N1=tensor(randn(sizet,sizet,sizet));
    N2=tensor(randn(sizet,sizet,sizet));
    nZ=norm(Z);
    nN1=norm(N1);

    % modify Z with the two different types of noise
    Zprime=Z+1/sqrt(100/l1-1)*nZ/nN1*N1;
    nZprime=norm(Zprime);

    N2Zprime=N2.*Zprime;
    nN2Zprime=norm(N2Zprime);

    Zdprime=Zprime+1/sqrt(100/l2-1)*nZprime/nN2Zprime*N2Zprime;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
