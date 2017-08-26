function out = cp_ngmres_o(Z,r,u0,par,par_ls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CP_NGMRES computes a CP factorization of a tensor using the N-GMRES
%optimization algorithm with ALS preconditioning
%
%   This function requires the Tensor toolbox (version 2.5) and the Poblano
%   toolbox. (http://www.sandia.gov/~tgkolda/TensorToolbox)
%
%   The nonlinear GMRES (N-GMRES) optimization method was proposed and analyzed
%   in the following references:
%      -Hans De Sterck, "A Nonlinear GMRES Optimization Algorithm for Canonical
%       Tensor Decomposition", SIAM J. Sci. Comput. 34, A1351-A1379, 2012.
%      -Hans De Sterck, "Steepest Descent Preconditioning for Nonlinear GMRES
%       Optimization", Numerical Linear Algebra with Applications 20, 453-471, 2013.
%
%   CP_NGMRES uses the ALS preconditioner as in the papers, and uses
%   poblano_linesearch from the Poblano toolbox. CP_NGMRES calls the N-GMRES
%   optimization method as implemented in ngmres.m.
%
%   out = CP_NGMRES(Z,r,u0,par,par_ls) computes an r-component CP factorization
%   of the tensor Z, via Frobenius norm minimization using the N-GMRES
%   optimization method.
%
%   INPUT:
%       Z: input tensor
%       r: number of components in the CP factorization
%       u0: cell array with the initial guess for the factor matrices
%       par: structure with the parameters for N-GMRES (see also ngmres.m)
%           par.w : maximum window size
%           par.maxIt : maximum number of iterations
%           par.relfTol : stopping tolerance on relative change in f
%           par.absgTol : stopping tolerance on the 2-norm of g
%           par.verbose : level of output to screen (0=no output, 1=summary output,
%               2=detailed output for each iteration)
%           par.logfev : flag to include the history of the number of f and g
%               evaluations in the output (true or false)
%           par.logf : flag to include the history of function values in the output
%               (true or false)
%           par.logg : flag to include the history of the 2-norm of the gradient
%               vectors in the output (true or false)
%           par.logRestart : flag to include the history of the restarts (true or
%               false)
%       par_ls: structure with parameters for poblano_linesearch (see also 
%               poblano_linesearch.m)
%           par_ls.LineSearch_ftol=1.000000000000000e-04;
%           par_ls.LineSearch_gtol=0.010000000000000;
%           par_ls.LineSearch_initialstep=1;
%           par_ls.LineSearch_maxfev=20;
%           par_ls.LineSearch_method='more-thuente';
%           par_ls.LineSearch_stpmax=1.000000000000000e+15;
%           par_ls.LineSearch_stpmin=1.000000000000000e-15;
%           par_ls.LineSearch_xtol=1.000000000000000e-15;
%
%   OUTPUT:
%       out: structure with result tensor and other output information
%           out.A : computed solution tensor
%           out.u: computed solution tensor in vector form
%           out.fev : history of number of f anf g evaluations (optional)
%           out.logf : history of function values (optional)
%           out.logg : history of the 2-norm of the gradient vectors (optional)
%           out.logRestart : restart history (indicator array) (optional)
%           out.convFlag : 1 = relfTol or absgTol convergence tolerances were
%                              satisfied
%                          2 = maxIt was reached without convergence
%                          0 = exit for other reason (NaN)
%
%   See also NGMRES, NGMRES_TEST_TENSOR_CP, TENSOR, SPTENSOR, KTENSOR, CP_OPT, POBLANO_LINESEARCH.
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for POBLANO
if ~exist('poblano_params','file')
    error(['CP_NGMRES requires Poblano Toolbox for Matlab. This can be ' ...
           'downloaded at http://software.sandia.gov/trac/poblano.']);
end

% Error checking
if ~isa(Z,'tensor')
    Z = tensor(Z);
end

if (nargin < 4)
    error('Error: invalid input arguments');
end

% CP factorization using CP_NGMRES
disp('+++ start n-gmres')
nZ2=norm(Z)^2;
out=ngmres_o(tt_fac_to_vec(u0), ...            % initial guess
           @(u) tt_cp_fun(u,Z,nZ2), ...      % function that computes f and g
           @(u,f,g) can_ALSu(Z,r,u,f,g,@(u) tt_cp_fun(u,Z,nZ2)), ...      
           ...                            % function that does one ALS iteration
           ...                            % (this is the preconditioner)
           @(fg,u0,f0,g0,d) poblano_linesearch(fg,u0,f0,g0,1.,d,par_ls),...
           ...                            % line search function
           par);                          % ngmres parameters
       
out.A=ktensor(tt_cp_vec_to_fac(out.u, Z));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,f,g,fev]=can_ALSu(T,r,u,f,g,fg);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function performs one ALS iteration (it uses a modification of the
%ALS routine in the tensor toolbox (only one iteration, and normalization
%is different))
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract number of dimensions and norm of T.
nModes = ndims(T);
dims=size(T);

% convert u to a ktensor
A = tt_cp_vec_to_fac(u,T);
A=ktensor(A);

A=can_ALS(T,r,A);

% convert A to a vector
u = tt_fac_to_vec(A.U);

% assume ALS is equivalent to 2 f and g evaluations (in reality it is a
% bit less expensive I think)
fev=2;

[f g]=fg(u);
fev=fev+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A=can_ALS(T,r,A);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function performs one ALS iteration
%this is a modification of the ALS routine in the tensor toolbox (only one
%iteration, and normalization is different)
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract number of dimensions and norm of T.
nModes = ndims(T);

% if there is a lambda factor in A, put it into the first mode
lambda=A.lambda;
A.U{1} = A.U{1} * spdiags(lambda,0,r,r);
A.lambda=ones(r,1);

% Iterate over all modes of the tensor
for n = 1:nModes
    % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
    Unew = mttkrp(T,A.U,n);

    % Compute the matrix of coefficients for linear system
    Y = ones(r,r);
    for i = [1:n-1,n+1:nModes]
        Y = Y .* (A.U{i}'*A.U{i});
    end

    Unew = (Y \ Unew')'; %<- Line from TTB 2.2.

    if issparse(Unew)
      A.U{n} = full(Unew);   % for the case R=1
    else
      A.U{n} = Unew;
    end
end

A = ktensor(A.U); % default: lambda=1

A=arrange(A); % normalize and introduce lambda, and make sure the components appear in order of lambda

% now distribute lambda evenly over all modes
facMat=spdiags(nthroot(A.lambda,nModes),0,r,r);
for n=1:nModes
    A.U{n} = A.U{n} * facMat;
end

A = ktensor(A.U); % lambda=1 by default
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
