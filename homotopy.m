function [ x, runtime ] = homotopy( x0, fun, N, hoption)
%This function implements the Homotopy method.  Choose x^(0) in R^n, and
%obtain the IVP system.
%Let lambda_i = i/N for i = 0, 1, …, N, solve the IVP system numerically at
%lambda_i using some IVP solvers, e.g., midpoint method and RK4.
%The stopping criterion for iterative methods is:
%Note that the backslash or the left matrix division operator is considered 
%in MATLAB to compute J(x)^-1*F(x).
%Inputs:
%x0 is the initial guess for the solution by default
%fun is a structure array by default containing the following fields:
%- .F is the given function F(x): R^n to R^n
%- .J is the Jacobian matrix J(x): R^n to R^nxn if available
%N is the last value of the iterations, starting at 0
%hoption specifies the IVP solver used in the homotopy method
%Outputs:
%x is the obtained approximate solution
%runtime is the running time for the method being used, which are obtained by using the MATLAB stopwatch timers tic and toc.
%Algorithm 10.4 Continuation Algorithm, to approximate the solution of the
%nonlinear system F(x)=0 given an initial approximation x

tic;
       
        yt=odesolver( 0, 1, x0, 1/N, fun, [], hoption);
        x=yt;

runtime=toc;

end