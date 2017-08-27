function [f g]=func_problemD(u)
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
    jodd = 1:2:n-1;
    jeven = 2:2:n;
    uTransformed = zeros(n,1);
    uTransformed(jodd) = 10*(u(jeven)-u(jodd).^2);
    uTransformed(jeven) = 1-u(jodd);

    % the value of the objective function
    f = 0.5*uTransformed'*uTransformed + 1;

    % the gradient
    g=zeros(n,1);
    g(jodd) = -20*u(jodd).*uTransformed(jodd) - uTransformed(jeven);
    g(jeven) = 10*uTransformed(jodd);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
