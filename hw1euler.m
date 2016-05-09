function [wc] = euler(ti, wi, h, f, argf)
% routine that implements Euler's method, algorithm 5.1 in the textbook
wc =wi + h*f(ti,wi,argf);
end