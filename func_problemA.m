function [f g]=func_problemA(u,d,u_exact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function implements optimization test problem A from "Steepest
%Descent Preconditioning for Nonlinear GMRES Optimization" by Hans De
%Sterck, arXiv:1106.4426v2, 2011
%
%OUTPUT: f is the value of the objective function, g is the gradient
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t=u-u_exact;

    % the value of the objective function
    f=0.5*t'*diag(d)*t + 1;

    % the gradient
    g=diag(d)*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
