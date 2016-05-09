function [fty] = fvdp(t, y, argf)
% this routine defines the f(t,y) in the Van der Pol problem
% the matrix size of the output argument fty is 1X2


dx = y(2);
dy = argf*(1 - y(1)^2)*y(2) - y(1);  %we plug in 1 and 2 because 

fty = [dx dy];
end