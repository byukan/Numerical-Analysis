function [ x, iter, runtime ] = steep( x0, fun, tol, maxiter )
%This function implements the Steepest Descent method.  We consider normalizing 
%the gradient to speed up.  See Algorithm 10.3 for the selection of alpha hat.
%Inputs:
%x0 is the initial guess for the solution by default
%fun is a structure array by default containing the following fields:
%- .F is the given function F(x): R^n to R^n
%- .J is the Jacobian matrix J(x): R^n to R^nxn if available
%tol is the tolerance
%maxiter is the maximal number of iterations allowed
%Outputs:
%x is the obtained approximate solution
%iter returns the number of iterations actually performed in those iterative methods 
%so that you are able to monitor the algorithm or even adjust maxiter.  It is defined 
%as ‘nan’ in the homotopy method
%runtime is the running time for the method being used, which are obtained by using the 
%MATLAB stopwatch timers tic and toc.
%Algorithm 10.3 Steepest Descent, to approximate a solution p to the
%minimization problem g(p)=min x in R^n of g(x) given an initial
%approximation x0

tic;
k=1;  %set up while counter step 1
x=x0;
while k<=maxiter  %step 2, do steps 3-15
    
    g=@(x) sum((fun.F(x)).^2);

	  %set approximate solution to initial guess as default

    g1=g(x);  %g1=g(x(k))
    
    Fx=fun.F(x);
    Jx=fun.J(x);

    z=2.*Jx'*Fx;  %z=gradient of (x(k))
    z0=norm(z,2);

    if z0==0
        x=[x;g1];  %zero gradient, the procedure completed, might have a minimum
        break;
    end
    
    z=z/z0;  %make z a unit vector
    a1=0;
    a3=1;
    g3=g(x-a3.*z);
    
    while g3>=g1  %step 6
        a3=a3/2;
        g3=g(x-a3.*z);
        if a3<tol/2
        x=[x;g(1)];  %no likely improvement, the procedure completed, might have a minimum
        break;
        end
    end
    
    a2=a3/2;  %step 9
    g2=g(x-a2.*z);
    
    h1=(g2-g1)/a2;  %step 10
    h2=(g3-g2)/(a3-a2);
    h3=(h2-h1)/a3;  %Newton's forward divided-difference formula is used to find the quadratic
    %P(a)=g1+h1*a+h3*a*(a-a2) that interpolates h(a) at a=0, a=a2, a=a3
    
    a0=.5*((a2-h1)/h3);  %the critical point of P occurs at a0
    g0=g(x-a0.*z);
    
    if g(x-a0.*z)==min(g0,g3)  %find a from a0, a3 so taht g=g(x-az)=min(g0,g3)
        a=a0;  %pick a0
    else if  g(x-a3.*z)==min(g0,g3)
         a=a3;  %pick a3
        end
    end
    g=min(g0,g3);
    
    x=x-a.*z;
    if abs(g-g1)<tol
        x=[x;g];  %the procedure was successful
    end
    iter=k;

    k=k+1;
end  %maximum number of iterations exceeded

runtime=toc;
    
end