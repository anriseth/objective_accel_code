function [out] = run_opts_tensor_CP(par)
    s = RandStream('mt19937ar','Seed', par.initSeed);
    RandStream.setGlobalStream(s);
  
    lineSearch_ngmres=@(fg,u0,f0,g0,d) poblano_linesearch(fg,u0,f0,g0,1.,d,par.par_ngmres);

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
    %% TODO: Save pars in par (from get_par)
    T=can_createTensorDense(pars);


    %--------------------------------------
    % do some pre-processing
    %--------------------------------------
    nT2=norm(T)^2;

    fg = @(u) tt_cp_fun(u,T,nT2);

    % Steepest descent preconditioner without linesearch
    M_sd=@(u,f,g) descentDir(u,f,g, ...
                             fg, ...                         % function that computes f and g
                             par.precStep1,par.precStep2);

    % ALS preconditioner
    M_als=@(u,f,g) can_ALSu(T,r,u,f,g,fg);

    %--------------------------------------
    % generate the random initial guess
    %--------------------------------------
    u0=rand(r*sum(size(T)),1);

    %--------------------------------------
    % call ngmres-als
    %--------------------------------------
    if par.compareNGMRESO_ALS==1
        out.out_ngmres_als=ngmres(u0, ...               % initial guess
                           fg, ...               % function that computes f and g
                           M_als, ...            % preconditioner function
                           lineSearch_ngmres,... % line search function
                           par.par_ngmres);      % ngmres parameters
    end

    %--------------------------------------
    % call other solvers for comparison
    %--------------------------------------
    if par.compareNGMRESO_ALS==1
        out.out_ngmreso_als=ngmres_o(u0, ...               % initial guess
                             fg, ...               % function that computes f and g
                             M_als, ...            % preconditioner function
                             lineSearch_ngmres,... % line search function
                             par.par_ngmres);      % ngmres parameters
    end
    if par.compareNGMRES_sd==1; % then N-GMRES, preconditioner: SD with small step
        disp('+++ start n-gmres with descent preconditioner')
        out.out_ngmres_sd=ngmres(u0, ...               % initial guess
                                 fg, ...               % function that computes f and g
                                 M_sd, ...             % preconditioner function
                                 lineSearch_ngmres,... % line search function
                                 par.par_ngmres);      % ngmres parameters
    end
    if par.compareNGMRESO_sd==1; % then N-GMRES-O, preconditioner: SD with small step
        disp('+++ start n-gmreso with descent preconditioner')
        out.out_ngmreso_sd=ngmres_o(u0, ...               % initial guess
                                    fg, ...               % function that computes f and g
                                    M_sd, ...             % preconditioner function
                                    lineSearch_ngmres,... % line search function
                                    par.par_ngmres);      % ngmres parameters
    end
    if par.compareNCG==1
        disp('+++ start n-cg directly')
        out.out_ncg=ncg(fg,u0,par.par_ncg);
        for i=2:size(out.out_ncg.TraceFuncEvals,2)
            out.out_ncg.TraceFuncEvals(i)=out.out_ncg.TraceFuncEvals(i)+out.out_ncg.TraceFuncEvals(i-1);
        end
    end
    if par.compareLBFGS==1
        disp('+++ start lbfgs')
        out.out_lbfgs=lbfgs(fg,u0,par.par_lbfgs);
        for i=2:size(out.out_lbfgs.TraceFuncEvals,2)
            out.out_lbfgs.TraceFuncEvals(i)=out.out_lbfgs.TraceFuncEvals(i)+out.out_lbfgs.TraceFuncEvals(i-1);
        end
    end
    if par.compareALS==1
        maxit = par.par_als.maxIt;
        fevALS = zeros(maxit,1);
        disp('+++ start ALS')
        u=u0;
        [f g]= fg(u);
        for k=1:maxit
            [u,f,g,fev]= M_als(u,f,g);%can_ALSu(T,r,u,f,g,fg);
            fALS(k)=f;
            gALS(k)=norm(g);
            fevALS(k) = fev;
        end
        for k = 2:maxit
            fevALS(k) = fevALS(k) + fevALS(k-1);
        end
        out.out_als.logf = fALS;
        out.out_als.logfev = fevALS;
    end
end
