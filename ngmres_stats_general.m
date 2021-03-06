function [outarr,iters,evals,fails] = ngmres_stats_general(numruns, ...
                                                      probnum,n,maxit)
    tol = 1e-10;
    outarr = cell(numruns,1);
    iters_ngmres_sdls = zeros(numruns,1);
    iters_ngmres_sd = zeros(numruns,1);
    iters_ngmreso_sdls = zeros(numruns,1);
    iters_ngmreso_sd = zeros(numruns,1);
    iters_ncg = zeros(numruns,1);
    iters_lbfgs = zeros(numruns,1);

    evals_ngmres_sdls = zeros(numruns,1);
    evals_ngmres_sd = zeros(numruns,1);
    evals_ngmreso_sdls = zeros(numruns,1);
    evals_ngmreso_sd = zeros(numruns,1);
    evals_ncg = zeros(numruns,1);
    evals_lbfgs = zeros(numruns,1);

    fails_ngmres_sdls = 0;
    fails_ngmres_sd = 0;
    fails_ngmreso_sdls = 0;
    fails_ngmreso_sd = 0;
    fails_ncg = 0;
    fails_lbfgs = 0;

    parfor seednum = 1:numruns
        fprintf('Iteration %d\n', seednum)
        par = get_par(seednum,probnum,n,maxit);
        outarr{seednum} = run_opts(par);

        fmin = par.probPars{probnum}(end);
        if isnan(fmin)
            fmin = min(outarr{seednum}.out_ngmres_sdls.logf);
            fmin = min(fmin, ...
                       min(outarr{seednum}.out_ngmres_sd.logf));
            fmin = min(fmin, ...
                       min(outarr{seednum}.out_ngmreso_sdls.logf));
            fmin = min(fmin, ...
                       min(outarr{seednum}.out_ngmreso_sd.logf));
            fmin = min(fmin, ...
                       min(outarr{seednum}.out_ncg.TraceFunc(2:end)));
            fmin = min(fmin, ...
                       min(outarr{seednum}.out_lbfgs.TraceFunc(2: ...
                                                              end)));
        end

        % Store the objective value at the initial condition
        f0 = outarr{seednum}.out_ncg.TraceFunc(1);

        % Store number of iterations to reach tolerance
        ind = find(outarr{seednum}.out_ngmres_sdls.logf < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_ngmres_sdls(seednum) = par.par_ngmres.maxIt;
            fails_ngmres_sdls = fails_ngmres_sdls + 1;
        else
            iters_ngmres_sdls(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ngmres_sd.logf < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_ngmres_sd(seednum) = par.par_ngmres.maxIt;
            fails_ngmres_sd = fails_ngmres_sd + 1;
        else
            iters_ngmres_sd(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ncg.TraceFunc(2:end) < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_ncg(seednum) = par.par_ngmres.maxIt;
            fails_ncg = fails_ncg + 1;
        else
            iters_ncg(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_lbfgs.TraceFunc(2:end) < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_lbfgs(seednum) = par.par_ngmres.maxIt;
            fails_lbfgs = fails_lbfgs + 1;
        else
            iters_lbfgs(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ngmreso_sdls.logf < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_ngmreso_sdls(seednum) = par.par_ngmres.maxIt;
            fails_ngmreso_sdls = fails_ngmreso_sdls + 1;
        else
            iters_ngmreso_sdls(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ngmreso_sd.logf < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_ngmreso_sd(seednum) = par.par_ngmres.maxIt;
            fails_ngmreso_sd = fails_ngmreso_sd + 1;
        else
            iters_ngmreso_sd(seednum) = ind;
        end

        % Store f/g evaluations to reach tolerance
        idx = min(length(outarr{seednum}.out_ngmres_sdls.logfev), ...
                  min(iters_ngmres_sdls(seednum)));
        evals_ngmres_sdls(seednum) = ...
            outarr{seednum}.out_ngmres_sdls.logfev(idx);
        idx = min(length(outarr{seednum}.out_ngmres_sd.logfev), ...
                  min(iters_ngmres_sd(seednum)));
        evals_ngmres_sd(seednum) = ...
            outarr{seednum}.out_ngmres_sd.logfev(idx);
        idx = min(length(outarr{seednum}.out_ncg.TraceFuncEvals), ...
                  1+iters_ncg(seednum));
        evals_ncg(seednum) = ...
            outarr{seednum}.out_ncg.TraceFuncEvals(idx);
        idx = min(length(outarr{seednum}.out_lbfgs.TraceFuncEvals), ...
                  1+iters_lbfgs(seednum));
        evals_lbfgs(seednum) = ...
            outarr{seednum}.out_lbfgs.TraceFuncEvals(idx);
        idx = min(length(outarr{seednum}.out_ngmreso_sdls.logfev), ...
                  min(iters_ngmreso_sdls(seednum)));
        evals_ngmreso_sdls(seednum) = ...
            outarr{seednum}.out_ngmreso_sdls.logfev(idx);
        idx = min(length(outarr{seednum}.out_ngmreso_sd.logfev), ...
                  min(iters_ngmreso_sd(seednum)));
        evals_ngmreso_sd(seednum) = ...
            outarr{seednum}.out_ngmreso_sd.logfev(idx);
    end

    iters.ngmres_sdls = iters_ngmres_sdls;
    iters.ngmres_sd = iters_ngmres_sd;
    iters.ngmreso_sdls = iters_ngmreso_sdls;
    iters.ngmreso_sd = iters_ngmreso_sd;
    iters.ncg = iters_ncg;
    iters.lbfgs = iters_lbfgs;

    evals.ngmres_sdls = evals_ngmres_sdls;
    evals.ngmres_sd = evals_ngmres_sd;
    evals.ngmreso_sdls = evals_ngmreso_sdls;
    evals.ngmreso_sd = evals_ngmreso_sd;
    evals.ncg = evals_ncg;
    evals.lbfgs = evals_lbfgs;

    fails.ngmres_sdls = fails_ngmres_sdls;
    fails.ngmres_sd = fails_ngmres_sd;
    fails.ngmreso_sdls = fails_ngmreso_sdls;
    fails.ngmreso_sd = fails_ngmreso_sd;
    fails.ncg = fails_ncg;
    fails.lbfgs = fails_lbfgs;
end
