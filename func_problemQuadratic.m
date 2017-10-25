function [f g]=func_problemQuadratic(u,T,u_exact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard quadratic function
%OUTPUT: f is the value of the objective function, g is the gradient
%
%by Asbjørn Nilsen Riseth, October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t=u-u_exact;

    % the value of the objective function
    f=0.5*t'*T*t;

    % the gradient
    g=T*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
