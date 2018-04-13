function Zdprime=can_createTensorDense(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function creates a test tensor for CP decomposition
%(this is the test from Tomasi and Bro (2006) and Acar et al. (2011) for
%a ramdom dense tensor modified to get specified collinearity between the
%columns)
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sizet=params(1);
    c=params(2);
    rTrue=params(3);
    l1=params(4);
    l2=params(5);

    % first generate K
    K=ones(rTrue,rTrue)*c;
    for k=1:rTrue
        K(k,k)=1;
    end
    K;

    % get C as the cholesky factor of K
    C=chol(K);

    for n=1:3
        % generate a random matrix
        M=randn(sizet,rTrue);
        % ortho-normalize the columns of M, gives Q
        [Q,R] = qr(M,0);
        U{n}=Q*C;
    end

    Z=full(ktensor(U));

    % generate two random tensors
    N1=tensor(randn(sizet,sizet,sizet));
    N2=tensor(randn(sizet,sizet,sizet));
    nZ=norm(Z);
    nN1=norm(N1);

    % modify Z with the two different types of noise
    Zprime=Z+1/sqrt(100/l1-1)*nZ/nN1*N1;
    nZprime=norm(Zprime);

    N2Zprime=N2.*Zprime;
    nN2Zprime=norm(N2Zprime);

    Zdprime=Zprime+1/sqrt(100/l2-1)*nZprime/nN2Zprime*N2Zprime;
end