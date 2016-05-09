function [ yx ] = linshoot( a,b,ya,yb,h,p,q,r )
%implements Algorithm 11.1, linear shooting method to approximate the solution of the
%boundary-value problem:
%-y'' + p(x)y' + q(x)y + r(x) = 0, a<=x<=b, y(a)=alpha, y(b)=beta:
%note: equations 11.3 and 11.4 are written as first-order systems and,
%solved
%a and b are endpoints of BVP
%ya is boundary condition alpha
%yb is boundary condition beta
%h is the step size
%p,q,r are the respective function handles of p(x), q(x) and r(x) of BVP
%output yx is a vector storing the approximations of the solution at xi

fun.dfdy = p;  %rename to pass into odesolver
fun.dfy = q;
fun.r = r;

argf = [];
option = 'rk4linshoot';

yt = odesolver(a, b, ya, h, fun, argf, option);


n = (b-a)/h;  % step 1, n is number of subintervals

u=yt{1};  %rename
v=yt{2};
w1 = ya;  % step 5
z = (yb-u(1,n))/v(1,n);  %z is w(2,0)
x = a;

    for i = 1 : n  %step 6
        w1(i) = u(1,i)+z*v(1,i);  %w1 is the approximation to y(xi)
        x(i) = a+i*h;
    end

    yx=[ya;w1'];
end



%Test example on page 676 in chapter 11.1
%linshoot(1,2,1,2, .1, @(x)(-2/x), @(x)(2/x^2), @(x) sin(log(x))/x^2)


%Solve Exercise 11.1.4b
%linshoot(0,pi/4,0,0, pi/20, @(x)0, @(x)-4, @(x) cos(x))