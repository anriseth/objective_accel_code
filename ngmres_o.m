function out=ngmres_o(u0,fg,M,linesearch,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NGMRES Nonlinear GMRES optimization method
%%%%
% THIS IS THE REVISED ALGORITHM THAT MINIMIZES THE OBJECTIVE
% INSTEAD OF THE \ell_2-NORM OF THE GRADIENT
%
% CHANGES BY: ASBJØRN NILSEN RISETH, 2017
%             anriseth@gmail.com
%
%%%%
% out=NGMRES(u0,fg,M,linesearch,par) minimizes a function using the nonlinear GMRES
%   (N-GMRES) optimization method.
%
%   The method was proposed and analyzed in the following references:
%      -Hans De Sterck, "A Nonlinear GMRES Optimization Algorithm for Canonical
%       Tensor Decomposition", SIAM J. Sci. Comput. 34, A1351-A1379, 2012.
%      -Hans De Sterck, "Steepest Descent Preconditioning for Nonlinear GMRES
%       Optimization", Numerical Linear Algebra with Applications 20, 453-471, 2013.
%
% INPUT:
% u0: initial guess (column vector)
% [f,g] = fg(u): fg is a function handle to a function that computes
%   function value f and gradient vector g in point u
% [u_new,f_new,g_new,fev] = M(u_old,f,g): M is a function handle to a
%   function that computes a new approximation u_new from an old
%   approximation u_old; M is the preconditioner of the n-gmres method;
%   f and g are the function and gradient in u_old, and f_new and g_new in
%   u_new; fev is the number of f and g evaluations
% [u,f,g,step,fev] = linesearch(fg,u0,f0,g0,d): linesearch is a function handle to
%   a function that does a line search minimizing f starting from u0 in
%   search direction d; f, g and step are the function value, gradient
%   vector and step length at the computed point u; f0 and g0 are the
%   function and gradient vector at u0, and fev is the number of f and g
%   evaluations during the line search
% par: parameters for ngmres
%   par.w : maximum window size
%   par.maxIt : maximum number of iterations
%   par.relfTol : stopping tolerance on relative change in f
%   par.absgTol : stopping tolerance on the 2-norm of g ( scaled by 1/length(u))
%   par.verbose : level of output to screen (0=no output, 1=summary output,
%       2=detailed output for each iteration)
%   par.logfev : flag to include the history of the number of f and g
%       evaluations in the output (true or false)
%   par.logf : flag to include the history of function values in the output
%       (true or false)
%   par.logg : flag to include the history of the 2-norm of the gradient
%       vectors in the output (true or false)
%   par.logRestart : flag to include the history of the restarts (true or
%       false)
%
% OUTPUT:
% out: ngmres output
%   out.u : computed solution point
%   out.fev : history of number of f anf g evaluations (optional)
%   out.logf : history of function values (optional)
%   out.logg : history of the 2-norm of the gradient vectors (optional)
%   out.logRestart : restart history (indicator array) (optional)
%   out.convFlag : 1 = relfTol or absgTol convergence tolerances were
%                       satisfied
%                  2 = maxIt was reached without convergence
%                  0 = exit for other reason
%
%  See also NGMRES_TEST_GENERAL, CP_NGMRES, NGMRES_TEST_TENSOR_CP
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some fixed parameters used in the algorithm
epsi=1e-12; % parameter used in regularizing the normal equation system

% some initialization
n=size(u0,1); % number of components in the u vectors
U=zeros(n,par.w); % solution vectors in the window
R=zeros(n,par.w); % gradient vectors in the window
logRestart=zeros(par.maxIt,1);

cnt=0; % counter for the number of iterations
nfev=0; % counter for the number of f and g evaluations

u=u0; % initialize the current iterate, u
[f,g]=fg(u); % compute function value and gradient in u
nfev=nfev+1;
ng=norm(g);
numu = length(u0);

% check convergence criteria
finishedIt= cnt >= par.maxIt;
finishedTol= ng <= par.absgTol*numu ;
finishedCrash=false;

while ~(finishedIt || finishedTol || finishedCrash)
    % start a new window
    logRestart(cnt+1)=1;
    if par.verbose==2
        disp('*** restart')
    end
    U(:,1)=u; % the current iterate is the first vector in the new window
    R(:,1)=g; % also store its gradient
    Q(1,1)=u'*g;
    %Q(1,1)=g'*g;


    restart=0; % no restart for now
    curw=1; % current window size (starts at 1, grows up to par.w)

    % start a sequence of accelerated updates, build towards window size w, but exit if a restart is required
    cntLoc=0; % counter of number of iterations since the last restart
    while restart==0
        cnt=cnt+1; % the very first new aproximation has cnt=1 (the initial guess is not counted) (but it is part of the window for acceleration!)
        cntLoc=cntLoc+1;

        f_previous=f; % save the function value of the previous iterate for convergence test

        % STEP I: get a new unaccelerated iterate
        %-------------------------------------------------
        [u_new,f_new,g_new,fev]=M(u,f,g);
        nfev=nfev+fev;

        % TODO: For practical implementations, include tolerance
        % check here (doesn't matter for the numerical experiments)

        % STEP II: compute the n-gmres accelerated iterate
        %-------------------------------------------------
        % compute its function value and gradient vector
        % form the least-squares system
        % efficiently by storing some previous inner products; see Washio
        % and Oosterlee, ETNA 6, pp. 271-290, 1997

        %eta=g_new'*g_new;
        eta=u_new'*g_new;
        for i=1:curw % iterations in window
            ksi1(i) = g_new'*U(:,i);
            ksi2(i) = u_new'*R(:,i);
            %ksi(i)=g_new'*R(:,i);
            beta(i,1)=eta-ksi1(i);
            %beta(i,1)=eta-ksi(i);
        end
        for i=1:curw % iterations in window
            for j=1:curw % iterations in window
                Mat(i,j)=Q(i,j)-ksi1(i)-ksi2(j)+eta;
            end
        end

        %delta=epsi*max(diag(Mat)+1e-16);
        delta=epsi*max(max(diag(Mat)),epsi); % the outer max is to avoid delta=0, which may occur if Mat=0, e.g. at numerical convergence

        % solve the normal equation system
        alpha=(Mat(1:curw,1:curw)+delta*eye(size(Mat(1:curw,1:curw),1))) \ beta(1:curw);

        % we could alternatively use Matlab's \ to solve the LS system...
        %matlabLS=0;
        %if matlabLS==1
        %    alpha=-Z\g_new;
        %end

        if isnan(norm(alpha)) % this does not seem to happen, but just in case
            finishedCrash=true;
            disp('+++++ WARNING: ngmres exit due to NaN')
        end

        % compute the accelerated approximation
        coef=1-sum(alpha);
        u_a=coef*u_new;

        for w=1:curw
            u_a=u_a+alpha(w)*U(:,w);
        end

        % STEP III: do a line search for globalization
        %-------------------------------------------------
        d=u_a-u_new; % search directory for the line search
                     % we'll do a restart if d is not a descent direction:
                     % (d is a descent direction if d * g_new < 0)
                     % (recal that g_new is the gradient, which is the direction of steepest ascent)
        if d' * g_new >= 0
            restart=1;
            if isempty(linesearch)
                u = u_new;
                f = f_new;
                g = g_new;
                step = 0.0;
                fev = 0;
            else
                [u,f,g,step,fev] = linesearch(fg,u_new,f_new,g_new, ...
                                              d);
            end
        else
            u = u_a;
            [f,g] = fg(u);
            step = 1.0;
            fev = 1;
        end

        nfev=nfev+fev;

        % provide some output and get some log information
        ng=norm(g);
        if par.verbose==2
            fprintf('it %d g %16.16g fev %d f %16.16g step %16.16g alpha %16.16g g_new %16.16g f_new %16.16g \n',...
                    cnt,ng,fev,f,step,norm(alpha),norm(g_new),f_new);
        end

        if par.logfev
            logfev(cnt)=nfev;
        end

        if par.logf
            logf(cnt)=f;
        end

        if par.logg
            logg(cnt)=ng;
        end

        if curw<par.w
            curw=curw+1;
        end

        if restart==0 % no restart, so update U and R, and Q (to store previous inner products)
                      % which column j of U and R to update
            j=mod(cntLoc,par.w)+1;
            U(:,j)=u;
            R(:,j)=g;
            for i=1:curw
                Q(i,j) = U(:,i)'*g;
                Q(j,i) = u'*R(:,i);
                %Q(j,i)=g'*R(:,i);
                %Q(i,j)=Q(j,i);
            end
        end

        % check the stopping criteria
        finishedIt= cnt >= par.maxIt;
        finishedTol= (ng <= par.absgTol*numu) || (abs(f-f_previous)/f_previous <= par.relfTol) ;
        if finishedIt || finishedTol || finishedCrash
            break % break out of the inner while loop when one of the convergence criteria is satisfied
        end

    end
end

% prepare the output data
out.u=u;
if par.logfev
    out.logfev=logfev;
end
if par.logf
    out.logf=logf;
end
if par.logg
    out.logg=logg;
end
if par.logRestart==true
    logRestart=logRestart(1:cnt);
    out.logRestart=logRestart;
end

out.convFlag=0; % finished because of crash
if finishedTol
    out.convFlag=1; % one of the convergence tolerances was reached
elseif finishedIt
    out.convFlag=2; % finished because maximum number of iterations was reached
end

if par.verbose==1
    fprintf('ngmres finished: it %d g %16.16g f %16.16g\n',cnt,ng,f);
end
