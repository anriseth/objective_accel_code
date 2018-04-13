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