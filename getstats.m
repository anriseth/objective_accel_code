function [outarr,iters,evals] = getstats(prob)
    probnum = prob.probnum;
    n = prob.n;
    maxits = prob.maxits;
    probname = prob.probnames{probnum};
    numruns = prob.numruns;
    qlevels = prob.qlevels;
    names = prob.names;
    runopt = prob.runopt;
    order = prob.order;

    fname = sprintf('data/problem%s%d',probname,n);
    if runopt == true
        [outarr,iters,evals] = ngmres_stats_general(numruns,probnum,n,maxits);
        save(strcat(fname,'.mat'), 'outarr', 'iters', 'evals');
    else
        load(strcat(fname,'.mat'));
    end

    qarr = quantlevels(evals, qlevels);
    tstr = latextable(qarr,names,order);
    fid = fopen(strcat(fname,'.tex'), 'w');
    fprintf(fid, tstr);
    fclose(fid);
end
