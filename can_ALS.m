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