function [f g]=func_problemG(u)
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
    
    t = sqrt(1e-5)*(u-1);
    t(n+1) = sum(u.^2)-0.25;
    % TODO: find out what the minimum value is
    f = 0.5*t'*t + 1; 
    
    
    % the gradient
    %g = sum(t)*sin(u) - t.*((1:n)'.*sin(u)+cos(u));
    g = sqrt(1e-5)*t(1:end-1) + 2*u.*t(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
