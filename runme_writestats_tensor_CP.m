runtensor_50 = true;
prob.names = {'N-GMRES-ALS', 'N-CG', 'L-BFGS', 'O-ACCEL-ALS', 'ALS'};
% We permute it as follows:
prob.order = [5, 4, 2, 3, 1];

% Report the following quantiles:
prob.qlevels = [0.1,0.5,0.9];
prob.numruns = 1000;

prob.maxits=1500;
prob.runopt=runtensor_50;
if prob.runopt
    [outarr,iters,evals,fails] = getstats_tensor_CP(prob);
end
