function [f g]=func_problemF(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function implements optimization test problem F from "Steepest
%Descent Preconditioning for Nonlinear GMRES Optimization" by Hans De
%Sterck, arXiv:1106.4426v2, 2011
%
%OUTPUT: f is the value of the objective function, g is the gradient
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = length(u);
    scosu = sum(cos(u));

    % NOTE: De Sterck uses a minus sign in front of (1:n),
    %       however, the original function in Moré et. al
    %       uses a plus sign.
    %       We use the plus sign
    t = n-scosu+(1:n)'.*(1-cos(u))-sin(u);
    f = 0.5*t'*t;

    % the gradient
    g = sum(t)*sin(u) + t.*((1:n)'.*sin(u)-cos(u));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
