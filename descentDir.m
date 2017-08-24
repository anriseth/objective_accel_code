function [u,f,g,fev]=descentDir(u,f,g,fg,precStep1,precStep2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function performs one steepest descent step with a small step length
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do a small step in the direction of steepest descent
    ng=norm(g);
    u=u-min(precStep1,precStep2*ng)*g/ng;

    [f g]=fg(u);
    fev=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
