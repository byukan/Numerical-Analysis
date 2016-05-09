function [ p ] = newton( p0, f, df, tol, itermax)
%if we are able to find the explicit expression of F', consider using
%Newton's method, Algorithm 2.3, for faster convergence

%this function also determines if a root is found, and the number of
%iterations actually used

i=1;
while i<=itermax
  p=p0-f(p0)/df(p0);
  if abs(p-p0)<tol
        found=true;
        p0=p;  %the procedure was successful
        i=itermax;
    end
    i=i+1;
    p0=p;  %update p0
end
p0=p;
found=false;  %the method failed
iter=i;  %iter is the number of iterations actually used
end