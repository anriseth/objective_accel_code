function [u,f,g,fev]=sd(u,f,g,fg,lsearch,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function performs one steepest descent step with line search
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do one steepest descent iteration
    d=-g;
    % De Sterck originally had a bug here.
    % He used
    % [u,f,g,step,fev] = lsearch(fg,u,f,g,d);
    % The fifth output from the poblano linesearch, however,
    % is the info number from More Thuente.
    % Which often returns 1, even though fev > 1
    [u,f,g,step,~,fev] = lsearch(fg,u,f,g,d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
