function [ x, iter, runtime ] = fsolver( x0, fun, paras )
%This function implements the solver functions in this program.  It takes 
%in arguments that make up the nonlinear system F(x)=0, then passes them to
%the actual solver functions.  Let F, a function from R^n to R^n, be 
%continuously differentiable with Jacobian matrix J(x).  This function solves
%the nonlinear system F(x)=0, using Newton’s method, Broyden’s method, Steepest
%Descent method, and Homotopy method.
%Inputs:
%x0 is the initial guess for the solution by default
%fun is a structure array by default containing the following fields:
%- .F is the given function F(x): R^n to R^n
%- .J is the Jacobian matrix J(x): R^n to R^nxn if available
%paras is a structure array containing the following fields:
%- .option specifies the method being used, which can be strings ‘newton’, ‘broyden’, ‘steep’, and ‘homotopy’
%- .tol is the tolerance if available
%- .maxiter is the maximal number of iterations allowed if available
%- .N is only used in the homotopy method, and empty in other methods
%- .hoption specifies the IVP solver used in the homotopy method, and empty in other methods
%
%Outputs:
%x is the obtained approximate solution
%iter returns the number of iterations actually performed in those iterative methods 
%so that you are able to monitor the algorithm or even adjust maxiter.  It is defined 
%as ‘nan’ in the homotopy method
%runtime is the running time for the method being used, which are obtained by using 
%the MATLAB stopwatch timers tic and toc.  To make the comparison fair, the following two lines:
%tic;
%runtime = toc;
%should be placed in fsolver at specific lines where each method starts and ends.  
%In addition, the stopping criterion should be placed somewhere specific in order 
%to get an approximate solution of the desired accuracy if maxiter is large enough.

switch paras.option
    case 'newton'
        
       [ newton_x, newton_iter, newton_runtime ] = newton( [0;0;0], fun, paras.tol, paras.maxiter )
        
    case 'broyden'
        
        [ broyden_x, broyden_iter, broyden_runtime ] = broyden( [0;0;0], fun, paras.tol, paras.maxiter )

    case 'steep'
        
        [ steep_x, steep_iter, steep_runtime ] = steep( [0;0;0], fun, paras.tol, paras.maxiter )

        
    case 'homotopy'
        
        [ homotopy_x, homotopy_runtime ] = homotopy( [0;0;0], fun, 10, 'rk4')

        
                

end