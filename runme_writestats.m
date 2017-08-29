%% 
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
run6_200    = true;
run6_500    = true;
run7_100    = false;
run7_200    = false;

prob.probnames{3} = 'B';
prob.probnames{4} = 'D';
prob.probnames{5} = 'E';
prob.probnames{6} = 'F';
prob.probnames{7} = 'G';

% All code runs things in this order:
prob.names = {'N-GMRES-ls', 'N-GMRES', 'NCG', 'LBFGS', 'N-GMRES-O-ls', 'N-GMRES-O'};
% We permute it as follows:
prob.order = [6, 5, 3, 4, 2, 1];

% Report the following quantiles:
prob.qlevels = [0.1,0.5,0.9];
prob.numruns = 1000;

%% 
% Problem B, n = 100
prob.maxits=1500;
prob.probnum=3;
prob.n = 100;
prob.runopt = run3_100;
getstats(prob);

%%
% Problem B, n = 200
prob.maxits=1500;
prob.probnum=3;
prob.n = 200;
prob.runopt = run3_200;
getstats(prob);
%%
% Problem D, n = 500
prob.maxits=500;
prob.probnum=4;
prob.n = 500;
prob.runopt = run4_500;
getstats(prob);
%%
% Problem D, n = 1000
prob.maxits=500;
prob.probnum=4;
prob.n = 1000;
prob.runopt = run4_1000;
getstats(prob);
%%
% Problem D, n = 50000
prob.maxits=500;
prob.probnum=4;
prob.n = 50000;
prob.runopt = run4_50000;
getstats(prob);
%%
% Problem D, n = 100000
prob.maxits=500;
prob.probnum=4;
prob.n = 100000;
prob.runopt = run4_100000;
getstats(prob);
%%
% Problem E, n = 100
prob.maxits=500;
prob.probnum=5;
prob.n = 100;
prob.runopt = run5_100;
getstats(prob);
%%
% Problem E, n = 200
prob.maxits=500;
prob.probnum=5;
prob.n = 200;
prob.runopt = run5_200;
getstats(prob);
%%
% Problem E, n = 50000
prob.maxits=500;
prob.probnum=5;
prob.n = 50000;
prob.runopt = run5_50000;
getstats(prob);
%%
% Problem E, n = 100000
prob.maxits=500;
prob.probnum=5;
prob.n = 100000;
prob.runopt = run5_100000;
getstats(prob);
%%
% Problem F, n = 200
prob.maxits=500;
prob.probnum=6;
prob.n = 200;
prob.runopt = run6_200;
getstats(prob);
%%
% Problem F, n = 500
prob.maxits=500;
prob.probnum=6;
prob.n = 500;
prob.runopt = run6_500;
getstats(prob);
%%
