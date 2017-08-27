function [outarr,iters,evals] = getstats(probnum,n,maxits,runopt)
    fname = sprintf('data/problem%s%d',probnames{probnum},n);
    if runopt == true
        [outarr,iters,evals] = ngmres_stats_general(numruns,probnum,n,maxits);
        save(strcat(fname,'.mat'), 'outarr', 'iters', 'evals');
    else
        load(fname);
    end

    qarr = quantlevels(evals, qlevels);
    tstr = latextable(qarr,names,order);
    fid = fopen(strcat(fname,'.tex'), 'w');
    fprintf(fid, tstr);
    fclose(fid);
end
