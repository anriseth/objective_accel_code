function [outarr,iters,evals,fails] = getstats_tensor_CP(prob)
    maxits = prob.maxits;
    numruns = prob.numruns;
    qlevels = prob.qlevels;
    names = prob.names;
    runopt = prob.runopt;
    order = prob.order;

    fprintf('*******\n Run tensor CP problem \n********\n');
    fname = sprintf('data/tensor_CP');
    if runopt == true
        [outarr,iters,evals,fails] = ngmres_stats_tensor_CP(numruns,maxits);
        save(strcat(fname,'.mat'), 'outarr', 'iters', 'evals','fails');
    else
        load(strcat(fname,'.mat'));
        if ~exist('fails')
            fails = 0;
        end
        if ~exist('outarr') % Outarr is too large for some problems
            outarr = NaN;
        end
    end

    qarr = quantlevels_tensor_CP(evals, qlevels);
    tstr = latextable(qarr,names,order);
    fid = fopen(strcat(fname,'.tex'), 'w');
    fprintf(fid, tstr);
    fclose(fid);
end
