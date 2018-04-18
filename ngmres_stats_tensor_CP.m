function [outarr,iters,evals,fails] = ngmres_stats_tensor_CP(numruns, maxit)
    tol = 1e-10;
    outarr = cell(numruns,1);
    iters_ngmres_als = zeros(numruns,1);
    iters_ngmreso_als = zeros(numruns,1);
    iters_als = zeros(numruns,1);
    iters_ncg = zeros(numruns,1);
    iters_lbfgs = zeros(numruns,1);

    evals_ngmres_als = zeros(numruns,1);
    evals_ngmreso_als = zeros(numruns,1);
    evals_als = zeros(numruns,1);
    evals_ncg = zeros(numruns,1);
    evals_lbfgs = zeros(numruns,1);

    fails_ngmres_als = 0;
    fails_ngmreso_als = 0;
    fails_als = 0;
    fails_ncg = 0;
    fails_lbfgs = 0;

    parfor seednum = 1:numruns
        fprintf('Iteration %d\n', seednum)
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
        %%

        outarr{seednum} = run_opts_tensor_CP(par);

        fmin = min(outarr{seednum}.out_ngmres_als.logf);
        fmin = min(fmin, ...
                   min(outarr{seednum}.out_ngmreso_als.logf));
        fmin = min(fmin, ...
                   min(outarr{seednum}.out_als.logf));
        fmin = min(fmin, ...
                   min(outarr{seednum}.out_ncg.TraceFunc(2:end)));
        fmin = min(fmin, ...
                   min(outarr{seednum}.out_lbfgs.TraceFunc(2:end)));

        % Store the objective value at the initial condition
        f0 = outarr{seednum}.out_ncg.TraceFunc(1);

        % Store number of iterations to reach tolerance
        ind = find(outarr{seednum}.out_ngmres_als.logf < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_ngmres_als(seednum) = par.par_ngmres.maxIt;
            fails_ngmres_als = fails_ngmres_als + 1;
        else
            iters_ngmres_als(seednum) = ind;
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
        ind = find(outarr{seednum}.out_ngmreso_als.logf < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_ngmreso_als(seednum) = par.par_ngmres.maxIt;
            fails_ngmreso_als = fails_ngmreso_als + 1;
        else
            iters_ngmreso_als(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_als.logf < (fmin+tol*(f0-fmin)),1);
        if isempty(ind)
            iters_als(seednum) = par.par_ngmres.maxIt;
            fails_als = fails_als + 1;
        else
            iters_als(seednum) = ind;
        end

        % Store f/g evaluations to reach tolerance
        idx = min(length(outarr{seednum}.out_ngmres_als.logfev), ...
                  min(iters_ngmres_als(seednum)));
        evals_ngmres_als(seednum) = ...
            outarr{seednum}.out_ngmres_als.logfev(idx);
        idx = min(length(outarr{seednum}.out_ncg.TraceFuncEvals), ...
                  1+iters_ncg(seednum));
        evals_ncg(seednum) = ...
            outarr{seednum}.out_ncg.TraceFuncEvals(idx);
        idx = min(length(outarr{seednum}.out_lbfgs.TraceFuncEvals), ...
                  1+iters_lbfgs(seednum));
        evals_lbfgs(seednum) = ...
            outarr{seednum}.out_lbfgs.TraceFuncEvals(idx);
        idx = min(length(outarr{seednum}.out_ngmreso_als.logfev), ...
                  min(iters_ngmreso_als(seednum)));
        evals_ngmreso_als(seednum) = ...
            outarr{seednum}.out_ngmreso_als.logfev(idx);
        idx = min(length(outarr{seednum}.out_als.logfev), ...
                  min(iters_als(seednum)));
        evals_als(seednum) = ...
            outarr{seednum}.out_als.logfev(idx);
    end

    iters.ngmres_als = iters_ngmres_als;
    iters.ngmreso_als = iters_ngmreso_als;
    iters.als = iters_als;
    iters.ncg = iters_ncg;
    iters.lbfgs = iters_lbfgs;

    evals.ngmres_als = evals_ngmres_als;
    evals.ngmreso_als = evals_ngmreso_als;
    evals.als = evals_als;
    evals.ncg = evals_ncg;
    evals.lbfgs = evals_lbfgs;

    fails.ngmres_als = fails_ngmres_als;
    fails.ngmreso_als = fails_ngmreso_als;
    fails.als = fails_als;
    fails.ncg = fails_ncg;
    fails.lbfgs = fails_lbfgs;
end
