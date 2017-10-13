function [f g]=func_problemE(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function implements optimization test problem E from "Steepest
%Descent Preconditioning for Nonlinear GMRES Optimization" by Hans De
%Sterck, arXiv:1106.4426v2, 2011
%
%OUTPUT: f is the value of the objective function, g is the gradient
%
%by Asbjørn Nilsen Riseth, September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = length(u);
    j1 = 1:4:n-3;
    j2 = 2:4:n-2;
    j3 = 3:4:n-1;
    j4 = 4:4:n;
    t = zeros(n,1);
    
    t(j1) = u(j1) + 10*u(j2);
    t(j2) = sqrt(5)*(u(j3)-u(j4));
    t(j3) = (u(j2) - 2*u(j3)).^2;
    t(j4) = sqrt(10)*(u(j1)-u(j4)).^2;
   
    % the value of the objective function
    f = 0.5*t'*t;

    % the gradient
    g=zeros(n,1);
    g(j1) = t(j1) + 2*sqrt(10)*(u(j1)-u(j4)).*t(j4);
    g(j2) = 10*t(j1) + 2*(u(j2)-2*u(j3)).*t(j3);
    g(j3) = sqrt(5)*t(j2) - 4*(u(j2)-2*u(j3)).*t(j3);
    g(j4) = -sqrt(5)*t(j2) - 2*sqrt(10)*(u(j1)-u(j4)).*t(j4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
