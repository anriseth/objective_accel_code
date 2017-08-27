function [out] = run_opts(par)
    s = RandStream('mt19937ar','Seed', par.initSeed);
    RandStream.setGlobalStream(s);

    %---------------------------------------------
    % set up the test problems
    %---------------------------------------------
    % the line search function for ngmres
    lineSearch_ngmres=@(fg,u0,f0,g0,d) poblano_linesearch(fg,u0,f0,g0,1.,d,par.par_ngmres);
    % the line search function for the sd preconditioner
    lineSearch_sdls=@(fg,u0,f0,g0,d) poblano_linesearch(fg,u0,f0,g0,1.,d,par.par_sdls);

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

    if par.problem==4
        % standard quadratic function with diagonal matrix A
        n=par.probPars{4}(1);
        u0=rand(n,1); % generate the initial guess
        fg=@(u) func_problemD(u); % set the objective function for N-GMRES
    end

    % the steepest descent preconditioner with line search
    M_sdls=@(u,f,g) sd(u,f,g, ...
                       fg, ...        % function that computes f and g
                       lineSearch_sdls,...
                       par.par_sdls); % sd parameters
                                      % the descent preconditioner with a fixed step
    M_sd=@(u,f,g) descentDir(u,f,g, ...
                             fg, ...                         % function that computes f and g
                             par.precStep1,par.precStep2);

    %--------------------------------------
    % call the ngmres method
    %--------------------------------------
    if par.compareNGMRES_sdls==1; % first N-GMRES, preconditioner: SD with line search
        disp('+++ start n-gmres with sdls preconditioner')
        out.out_ngmres_sdls=ngmres(u0, ...               % initial guess
                                   fg, ...               % function that computes f and g
                                   M_sdls, ...           % preconditioner function
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
    if par.compareNGMRESO_sdls==1; % first N-GMRES, preconditioner: SD with line search
        disp('+++ start n-gmres-o with sdls preconditioner')
        out.out_ngmreso_sdls=ngmres_o(u0, ...               % initial guess
                                      fg, ...               % function that computes f and g
                                      M_sdls, ...           % preconditioner function
                                      lineSearch_ngmres,... % line search function
                                      par.par_ngmres);      % ngmres parameters
    end
    if par.compareNGMRES_sd==1; % then N-GMRES, preconditioner: SD with small step
        disp('+++ start n-gmreso with descent preconditioner')
        out.out_ngmreso_sd=ngmres_o(u0, ...               % initial guess
                                    fg, ...               % function that computes f and g
                                    M_sd, ...             % preconditioner function
                                    lineSearch_ngmres,... % line search function
                                    par.par_ngmres);      % ngmres parameters
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
    if par.compare_sdls==1
        disp('+++ start sdls preconditioner only')
        u=u0;
        [f g]=fg(u);
        nfev=0;
        for k=1:par.par_ngmres.maxIt
            [u,f,g,fev]=M_sdls(u,f,g);
            nfev=nfev+fev;
            out.out_sdls.logf(k)=f;
            out.out_sdls.logg(k)=norm(g);
            out.out_sdls.logfev(k)=nfev;
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
