function [f g]=func_problemC(u,T,u_exact,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function implements optimization test problem B from "Steepest
%Descent Preconditioning for Nonlinear GMRES Optimization" by Hans De
%Sterck, arXiv:1106.4426v2, 2011
%
%OUTPUT: f is the value of the objective function, g is the gradient
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = length(u);
    u=u-u_exact;
    
    cond(T)
    
    t=u;
    t(2:end)=t(2:end)-alpha*t(1)^2;
    
    % the value of the objective function
    f=0.5*t'*T*t + 1;

    % the gradient
    g=T*t;
    g(1) = g(1) - 2*alpha*u(1)*sum(g(2:end));
    %g(2:end)=g(2:end)-alpha*u(1)^2*d(2:end);
    %g(1)=g(1)-2*alpha*u(1)*d(2:end)'*u(2:end)+2*alpha^2*u(1)^3*sum(d(2:end));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
