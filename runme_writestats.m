%%
run1_100    = true;
run1_200    = false;
run2_100    = false;
run2_200    = false;
run3_100    = false;
run3_200    = false;
run4_500    = false;
run4_1000   = false;
run4_50000  = false;
run4_100000 = false;
run5_100    = false;
run5_200    = false;
run5_50000  = false;
run5_100000 = false;
run6_200    = false;
run6_500    = false;
run7_100    = false;
run7_200    = false;

prob.probnames{1} = 'A';
prob.probnames{2} = 'B';
prob.probnames{3} = 'C';
prob.probnames{4} = 'D';
prob.probnames{5} = 'E';
prob.probnames{6} = 'F';
prob.probnames{7} = 'G';

% All code runs things in this order:
prob.names = {'N-GMRES-A', 'N-GMRES-B', 'N-CG', 'L-BFGS', 'O-ACCEL-A', 'O-ACCEL-B'};
% We permute it as follows:
prob.order = [6, 5, 3, 4, 2, 1];

% Report the following quantiles:
prob.qlevels = [0.1,0.5,0.9];
prob.numruns = 1000;

%%
% Problem A, n = 100
prob.maxits=1500;
prob.probnum=1;
prob.n = 100;
prob.runopt = run1_100;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem A, n = 200
prob.maxits=1500;
prob.probnum=1;
prob.n = 200;
prob.runopt = run1_200;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem B, n = 100
prob.maxits=1500;
prob.probnum=2;
prob.n = 100;
prob.runopt = run2_100;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end

%%
% Problem B, n = 200
prob.maxits=1500;
prob.probnum=2;
prob.n = 200;
prob.runopt = run2_200;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem C, n = 100
prob.maxits=1500;
prob.probnum=3;
prob.n = 100;
prob.runopt = run3_100;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end

%%
% Problem C, n = 200
prob.maxits=1500;
prob.probnum=3;
prob.n = 200;
prob.runopt = run3_200;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem D, n = 500
prob.maxits=1500;
prob.probnum=4;
prob.n = 500;
prob.runopt = run4_500;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem D, n = 1000
prob.maxits=1500;
prob.probnum=4;
prob.n = 1000;
prob.runopt = run4_1000;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem D, n = 50000
% Problem D, n = 100000
% Are at the bottom of the file
%%
% Problem E, n = 100
prob.maxits=1500;
prob.probnum=5;
prob.n = 100;
prob.runopt = run5_100;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem E, n = 200
prob.maxits=1500;
prob.probnum=5;
prob.n = 200;
prob.runopt = run5_200;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem D, n = 50000
% Problem D, n = 100000
% Are at the bottom of the file
%%
% Problem F, n = 200
prob.maxits=1500;
prob.probnum=6;
prob.n = 200;
prob.runopt = run6_200;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem F, n = 500
prob.maxits=1500;
prob.probnum=6;
prob.n = 500;
prob.runopt = run6_500;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem G, n = 100
prob.maxits=1500;
prob.probnum=7;
prob.n = 100;
prob.runopt = run7_100;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem G, n = 200
prob.maxits=1500;
prob.probnum=7;
prob.n = 200;
prob.runopt = run7_200;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end

%% These take much longer
% Problem D, n = 50000
prob.maxits=1500;
prob.probnum=4;
prob.n = 50000;
prob.runopt = run4_50000;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem D, n = 100000
prob.maxits=1500;
prob.probnum=4;
prob.n = 100000;
prob.runopt = run4_100000;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem E, n = 50000
prob.maxits=1500;
prob.probnum=5;
prob.n = 50000;
prob.runopt = run5_50000;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
%%
% Problem E, n = 100000
prob.maxits=1500;
prob.probnum=5;
prob.n = 100000;
prob.runopt = run5_100000;
if prob.runopt
    [outarr,iters,evals,fails] = getstats(prob);
end
