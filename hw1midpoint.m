function [wc] = midpoint(ti, wi, h, f, argf)
% routine that implements the Midpoint Method.
wc=wi+h*f(ti+(h/2),wi+(h/2)*f(ti,wi,argf));
end 