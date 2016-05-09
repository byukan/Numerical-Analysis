function [ yx, error, actual ] = bvpsolver( xinit, yinit, fun, h, option )
%solves the two-point second order BVP (1)
%xinit is a vector containing the endpoints of the interval ie a and b
%yinit is a vector containing the function values at the endpoints ie alpha
%and beta
%fun is a structure array containing the following fields:
%.f the function f(x,y,y') which can be defined as fxydy=f(x,y,dy,argf) in
%a separate function m-file where argf contains related parameters, or an
%anonymous function f=@(x,y,dy,argf)(...); you may skip the argf if it is
%not necessary
%.dfdy the function fy'(x,y,y') which is reduced to P(x) in the linear case
%.dfy the function fy(x,y,y') which is reduced to q(x) in the linear case
%.r the function r(x) only in the linear case.  Note that you may skip this
%field fo rhte nonlinear case
%For example, if the linear BVP has f(x,y,y')=x^2*y'+x*y+x in (1), then you
%may define in the main.m
%fun.f=@(x,y,dy)(x.^2.*dy+x.*y+x);
%fun.dfdy=@(x) (x.^2);
%fun.dfy=@(x)(x);
%fun.r=@(x)(x);
%or:  fun=struct('f', @(x,y,dy)(x.^2.*dy+x.*y+x), 'dfdy', @(x)(x.^2),
%'dfy', @(x)(x), 'r', @(x)(x));
%the structure array simplifies syntax and groups arguments of similar
%attributes, which is very convenient in coding.
%h is the step size:  note that xi=a+ih for i from 0 to N+1, and b-a=(N+1)h
%when using the finite-difference method, while b-a=NJh when using the
%shotting method.  We determine N nodes in shotting methods, but N-1
%interior nodes in finite-difference methods.
%option is a string specifying hte method to be used which has the
%following choices: 'linshoot' - linear shooting method; 'nlinshoot' -
%nonlinear shooting method; 'lindiff' - linear finite-difference method;
%'nlindiff' - nonlinear finite-difference method
%yx is a vector storing the approximations of the soln at xi
%error is the global error defined as max i of abs(y(xi)-wi)

%To compute "error" in bvpsolver.m, we include one more field in the
%input argument "fun", e.g., "fun.y" which contains the actual solution y(t).


switch option
    case 'linshoot'  %linear shooting method
        
        a = xinit(1);
        b = xinit(2);
        ya = yinit(1);
        yb = yinit(2);
        p=fun.dfdy;
        q=fun.dfy;
        r=fun.r;
       
        yx = linshoot( a,b,ya,yb,h,p,q,r );
        
        n = (b-a)/h; %now check error, n is number of subintervals
        for i = 1:n
            x(i) = a+(i-1)*h;
            
            %actual solution y(t) to compute error
            fun.y(i) = (-1/3)*cos(2*x(i))-(sqrt(2)/6)*sin(2*x(i))+(1/3)*cos(x(i)); %Exercise 11.1.4b
            
            error(i)=abs(fun.y(i)-yx(i));  %fun.y(i) is the actual solution, and yx(i) is the approximation
            error=max(error);
        end
        
    case 'nlinshoot'  %nonlinear shooting method
        
        a = xinit(1);
        b = xinit(2);
        ya = yinit(1);
        yb = yinit(2);
        
        tol=10^-4;  %declare this manually for the exercises since there's no way to
        %pass it from odesolver, the way the function arguments are in the assignment
        itermax=10;
        
        yx = nlinshoot( a, b, ya, yb, h, fun, itermax, tol );
        
        n = (b-a)/h; %now check error, n is number of subintervals
        for i = 1:n
            x(i) = a+(i-1)*h;
            
            %actual solution y(t) to compute error
            fun.y(i) = x(i)+(1/x(i)); %Exercise 11.2.4b
            
            error(i)=abs(fun.y(i)-yx(i));  %fun.y(i) is the actual solution, and yx(i) is the approximation
            error=max(error);
        end
        
        
    case 'lindiff'  %linear finite-difference method
        
        a = xinit(1);
        b = xinit(2);
        ya = yinit(1);
        yb = yinit(2);
        p=fun.dfdy;
        q=fun.dfy;
        r=fun.r;
       
        yx = lindiff( a,b,ya,yb,h,p,q,r );
        
        n = (b-a)/h; % step 1, n is number of subintervals
        for i = 1:n
            x(i) = a+(i-1)*h;
            
            %actual solution y(t) to compute error
            fun.y(i) = (-1/3)*cos(2*x(i))-(sqrt(2)/6)*sin(2*x(i))+(1/3)*cos(x(i)); %Exercise 11.1.4b
            
            error(i)=abs(fun.y(i)-yx(i));  %fun.y(i) is the actual solution, and yx(i) is the approximation
            error=max(error);
        end
        
    case'nlindiff'  %nonlinear finite-difference method
        
        a = xinit(1);
        b = xinit(2);
        ya = yinit(1);
        yb = yinit(2);
        
        tol=10^-4;  %declare this manually for the exercises since there's no way to
        %pass it from odesolver, the way the function arguments are in the assignment
        itermax=10;
        
        yx = nlindiff( a, b, ya, yb, h, fun, itermax, tol );
        
        n = (b-a)/h; %now check error, n is number of subintervals
        for i = 1:n
            x(i) = a+(i-1)*h;
            
            %actual solution y(t) to compute error
            fun.y(i) = x(i)+(1/x(i)); %Exercise 11.2.4b
            
            error(i)=abs(fun.y(i)-yx(i));  %fun.y(i) is the actual solution, and yx(i) is the approximation
            error=max(error);
        end



end
        actual=[fun.y';yb];
end