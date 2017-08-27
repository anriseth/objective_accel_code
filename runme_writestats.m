%% 
run3_100    = false;
run3_200    = false;
run4_500    = false;
run4_1000   = false;
run4_50000  = false;
run4_100000 = false;
run5_100    = true;
run5_200    = true;
run5_50000  = true;
run5_100000 = true;

probnames{3} = 'B';
probnames{4} = 'D';
probnames{4} = 'E';
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
getstats(probnum,n,maxits,run3_100);

%%
% Problem B, n = 200
maxits=1500;
probnum=3;
n = 200;
getstats(probnum,n,maxits,run3_200);
%%
% Problem D, n = 500
maxits=500;
probnum=4;
n = 500;
getstats(probnum,n,maxits,run4_500);
%%
% Problem D, n = 1000
maxits=500;
probnum=4;
n = 1000;
getstats(probnum,n,maxits,run4_1000);
%%
% Problem D, n = 50000
maxits=500;
probnum=4;
n = 50000;
getstats(probnum,n,maxits,run4_50000);
%%
% Problem D, n = 100000
maxits=500;
probnum=4;
n = 100000;
getstats(probnum,n,maxits,run4_100000);
%%
% Problem E, n = 100
maxits=500;
probnum=5;
n = 100;
getstats(probnum,n,maxits,run5_100);
%%
% Problem E, n = 200
maxits=500;
probnum=5;
n = 200;
getstats(probnum,n,maxits,run5_200);
%%
% Problem E, n = 50000
maxits=500;
probnum=5;
n = 50000;
getstats(probnum,n,maxits,run5_50000);
%%
% Problem E, n = 100000
maxits=500;
probnum=5;
n = 100000;
getstats(probnum,n,maxits,run5_100000);
%%
