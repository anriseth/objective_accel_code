function [outarr,iters,evals] = ngmres_stats_general(numruns,probnum,n,maxit)
    fmin = 1.0;
    tol = 1e-6;
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

    parfor seednum = 1:numruns
        fprintf('Iteration %d\n', seednum)
        par = get_par(seednum,probnum,n,maxit);
        outarr{seednum} = run_opts(par);

        % Store number of iterations to reach tolerance
        ind = find(outarr{seednum}.out_ngmres_sdls.logf < fmin+tol,1);
        if isempty(ind)
            iters_ngmres_sdls(seednum) = par.par_ngmres.maxIt;
        else
            iters_ngmres_sdls(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ngmres_sd.logf < fmin+tol,1);
        if isempty(ind)
            iters_ngmres_sd(seednum) = par.par_ngmres.maxIt;
        else
            iters_ngmres_sd(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ncg.TraceFunc(2:end) < fmin+tol,1);
        if isempty(ind)
            iters_ncg(seednum) = par.par_ngmres.maxIt;
        else
            iters_ncg(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_lbfgs.TraceFunc(2:end) < fmin+tol,1);
        if isempty(ind)
            iters_lbfgs(seednum) = par.par_ngmres.maxIt;
        else
            iters_lbfgs(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ngmreso_sdls.logf < fmin+tol,1);
        if isempty(ind)
            iters_ngmreso_sdls(seednum) = par.par_ngmres.maxIt;
        else
            iters_ngmreso_sdls(seednum) = ind;
        end
        ind = find(outarr{seednum}.out_ngmreso_sd.logf < fmin+tol,1);
        if isempty(ind)
            iters_ngmreso_sd(seednum) = par.par_ngmres.maxIt;
        else
            iters_ngmreso_sd(seednum) = ind;
        end

        % Store f/g evaluations to reach tolerance
        evals_ngmres_sdls(seednum) = ...
            outarr{seednum}.out_ngmres_sdls ...
            .logfev(iters_ngmres_sdls(seednum));
        evals_ngmres_sd(seednum) = ...
            outarr{seednum}.out_ngmres_sd ...
            .logfev(iters_ngmres_sd(seednum));
        evals_ncg(seednum) = ...
            outarr{seednum}.out_ncg.TraceFuncEvals(1+ ...
                                                   iters_ncg(seednum));
        evals_lbfgs(seednum) = ...
            outarr{seednum}.out_lbfgs.TraceFuncEvals(1+ ...
                                                     iters_lbfgs(seednum));        
        evals_ngmreso_sdls(seednum) = ...
            outarr{seednum}.out_ngmreso_sdls ...
            .logfev(iters_ngmreso_sdls(seednum));
        evals_ngmreso_sd(seednum) = ...
            outarr{seednum}.out_ngmreso_sd ...
            .logfev(iters_ngmreso_sd(seednum));
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
end
