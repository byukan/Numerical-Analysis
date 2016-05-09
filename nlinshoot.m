function [ yx ] = nlinshoot( a, b, ya, yb, h, fun, itermax, tol )
%implements algorithm 11.2, nonlinear shooting method where newton's method is used
%to update tk
% To approximate the solution of the nonlinear boundary-value problem
%    y'' = f(x,y,y'), a<=x<=b, y(a) = ya, y(b) = ya:
%a and b are the endpoints of BVP
%ya, and yb are boundary conditions alpha and beta
%h is the step size used to calulate n, which is the number of subintervals
%fun is the same as in bvpsolver.m
%itermax is the maximal number of iterations allowed
%tol is the tolerance in the stopping critera abs(y(xn+1, tk)-yb) <tol
% output yx is the vector of approximations w(1,i) to y(xi) for each i=0,1,...,n



option = 'rk4nlinshoot';



n = (b-a)/h;  % step 1, n is number of subintervals
k=1;
TK = (yb-ya)/(b-a);  %TK could also be input
w1 = zeros(1,n+1);
w2 = zeros(1,n+1);

while k <= itermax % step 2, do steps 3-10

    w1(1) = ya;  % step 3
    w2(1) = TK;

    
    argf{1}=w1;  %pass w1 and w2 into odesolver
    argf{2}=w2;
    yt = odesolver(a, b, ya, h, fun, argf, option);

    w1=yt{1};  %make yt a cell to store vectors u and v to get passed back to linshoot.m
    u1=yt{2};
    
        if abs(w1(n+1)-yb) <= tol   % step 7, test for accuracy

            for i = 1 : n+1  % step 8
                x(i) = a+i*h;
            end
        end

        TK = TK-(w1(n+1)-yb)/u1;  % step 10, Newton's method is used to compute TK
        k = k+1;
        
end

yx=w1';

end
%Test example on page 683 in chapter 11.2
%fun.f = @(x,y,dy) -(y./8).*dy+x.^3./4+4;
%fun.dfdy = @(x,y,dy) -(y./8);
%fun.dfy = @(x,y,dy) -(dy./8);
%fun.r = @(x,y,dy) x.^3./4+4;
%nlinshoot( 1, 3, 17, 43/3, .1, fun, 10, 10^-5 )







%Now run the nonlinear BVP using nlinshoot and nlindiff
%fun.f = @(x,y,dy) 2*y^3-6*y-2*x^3;
%fun.dfdy = @(x,y,dy) 0;  %partial of f with respect to y'
%fun.dfy = @(x,y,dy) 6*y^2-6;  %partial of f with respect to y

%Exercise 11.2.4b
%[ yx2, error2 ] = bvpsolver( [1,2], [2,5/2], fun, .1, 'nlinshoot' )