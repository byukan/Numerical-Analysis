function [ x, iter, runtime ] = broyden( x0, fun, tol, maxiter )
%This function implements the Broyden’s method.
%Inputs:
%x0 is the initial guess for the solution by default
%fun is a structure array by default containing the following fields:
%- .F is the given function F(x): R^n to R^n
%- .J is the Jacobian matrix J(x): R^n to R^nxn if available
%tol is the tolerance
%maxiter is the maximal number of iterations allowed
%Outputs:
%x is the obtained approximate solution
%iter returns the number of iterations actually performed in those iterative methods so that you are able to monitor the algorithm or even adjust maxiter.  It is defined as ‘nan’ in the homotopy method
%runtime is the running time for the method being used, which are obtained by using the MATLAB stopwatch timers tic and toc.
%Algorithm 10.2 Broyden, to approximate the solution of the nonlinear
%system F(x)=0 given an initial approximation x

tic;  %start timer
iter=0;
x=x0;  %set approximate solution to initial guess as default

F = fun.F(x0);
J = fun.J(x0);

A0=inv(J);
dx=-A0 * F;

x0 = x0+dx;


for i=2:maxiter
    
    w = fun.F(x0);  %save v, note: v=F(x0)
    dy = w - F;
    u = A0 * dy;
    v = dx' * A0;
    p = dx' * u;
    A0 = A0 + ( dx - u ) * v / p;
    dx = -A0 * w;
    x0 = x0 + dx;

iter=iter+1;
    
    if (max(abs(dx))<tol) 
       
          y = x0;
              x=y;
          break
 
    else
       F = w;
    end
    
    runtime=toc;

    
end
 
end