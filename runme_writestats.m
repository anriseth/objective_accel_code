%% 
run3_100 = false;
run3_200 = false;
run4_500 = false;
run4_1000 = false;
run4_50000 = true;
run4_100000 = true;

probnames{3} = 'B';
probnames{4} = 'D';
% All code runs things in this order:
names = {'N-GMRES-ls', 'N-GMRES', 'NCG', 'LBFGS', 'N-GMRES-O-ls', 'N-GMRES-O'};
% We permute it as follows:
order = [6, 5, 3, 4, 2, 1];

% Report the following quantiles:
qlevels = [0.1,0.5,0.9];
numruns = 1000;
%% 
% Problem B, n = 100
maxits=1500;
probnum=3;
n = 100;
fname = sprintf('data/problem%s%d',probnames{probnum},n);
if run3_100 == true
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
%%
% Problem B, n = 200
maxits=1500;
probnum=3;
n = 200;
fname = sprintf('data/problem%s%d',probnames{probnum},n);
if run3_200 == true
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
%%
% Problem D, n = 500
maxits=500;
probnum=4;
n = 500;
fname = sprintf('data/problem%s%d',probnames{probnum},n);
if run4_500 == true
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
%%
% Problem D, n = 1000
maxits=500;
probnum=4;
n = 1000;
fname = sprintf('data/problem%s%d',probnames{probnum},n);
if run4_1000 == true
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
%%
% Problem D, n = 50000
maxits=500;
probnum=4;
n = 50000;
fname = sprintf('data/problem%s%d',probnames{probnum},n);
if run4_50000 == true
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
%%
% Problem D, n = 100000
maxits=500;
probnum=4;
n = 100000;
fname = sprintf('data/problem%s%d',probnames{probnum},n);
if run4_50000 == true
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
%%