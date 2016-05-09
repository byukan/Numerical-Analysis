function [ yx ] = nlindiff( a, b, ya, yb, h, fun, itermax, tol )
%applies the finite-difference method and the Newton's method for nonlinear
%systems to solve the nonlinear BVP
%a and b are endpoints of BVP
%ya is boundary condition alpha
%yb is boundary condition beta
%itermax and tol are the maximal number of itermations allowed and the
%tolernace respectively in the Newton's method for nonliner systems
%the Jacobian matrix J cna be constructed similar to A in the linear case.
%Be careful that J varies at each iteration which involves the function fy,
%and fy'.  Similar to lindiff, instead of solving a linear system with ethe
%coefficient matrix J, the matrix left division '\' can be used.
%the initial values for wi for i from 1 to N are chosen as what are used
%in Algorithm 11.4 on page 693 step 2

%output yx is a vector storing the approximations of the solution at xi


n=(b-a)/h;  %difference methods number of subintervals, use same formula for this algorithm



w=zeros(n-1,1);

for i=1:n-1;  %create vector w of the iterates that converge to the solution of system
    w(i)= ya+i*((yb-ya)/(b-a))*h;
end

k=1;  %step 3

while k<=itermax
  
%create the first entries of d1, d2, and d3, since t has to evaluate at w(i-1) in the next for loop
x(1)=a+h;  %x is the vector of time steps used to calculate entries of d1, d2, and d3
        t(1)=(w(2)-ya)/(2*h);
        d1(1)=-1-(h/2).*fun.dfdy(x(1),w(1),t(1));  %d1 is the lower diagonal of tridiagonal matrix
    	d3(1)=-1+(h/2).*fun.dfdy(x(1),w(1),t(1));  %d3 is the upper diagonal of tridiagonal matrix
        d2(1)=2+h^2.* fun.dfy(x(1),w(1),t(1));  %d2 is the main diagonal of tridiagonal matrix
    
    for i=2:n-2  %sparse, spdiags will automatically cut off d1(1) and d3(n), that's how it
    %works, so make all the diagonals the same size.
        x(i)=a+i*h;  %x is the vector of time steps used to calculate entries of d1, d2, and d3
        t(i)=(w(i+1)-w(i-1))/(2*h);
        d1(i)=-1-(h/2).*fun.dfdy(x(i),w(i),t(i));  %d1 is the lower diagonal of tridiagonal matrix
    	d3(i)=-1+(h/2).*fun.dfdy(x(i),w(i),t(i));  %d3 is the upper diagonal of tridiagonal matrix
        d2(i)=2+h^2.* fun.dfy(x(i),w(i),t(i));  %d2 is the main diagonal of tridiagonal matrix
    end
    
    
x(n-1)=a+(n-1)*h;  %x is the vector of time steps used to calculate entries of d1, d2, and d3
t(n-1)=(yb-w(n-2))/(2*h);
d1(n-1)=-1-(h/2).*fun.dfdy(x(n-1),w(n-1),t(n-1));  %d1 is the lower diagonal of tridiagonal matrix
d3(n-1)=-1+(h/2).*fun.dfdy(x(n-1),w(n-1),t(n-1));  %d3 is the upper diagonal of tridiagonal matrix
d2(n-1)=2+h^2.* fun.dfy(x(n-1),w(n-1),t(n-1));  %d2 is the main diagonal of tridiagonal matrix
    
    
b(1)= -(2*w(1)-w(2)-ya+h^2*fun.f(x(1),w(1),t(1)));  %%construct vector b, starting with first entry
b(n-1)= -(2*w(n-1)-w(n-2)-yb+h^2*fun.f(x(n-1),w(n-1),t(n-1)));  %last entry
    for i=2:n-2  %fill in inner entries of b
        b(i)= -(2*w(i)-w(i+1)-w(i-1)+h^2*fun.f(x(i),w(i),t(i)));
    end

J=spdiags([d3',d2',d1'], -1:1, n-1, n-1);  %J is the Jacobian matrix
%to make spdiags work, I had to take the transpose so the size works with
%this function, also this function is written to use n-1 as size

if norm(J'\b',2)< tol
    break;
end

k=k+1;
w = w+J'\b';
end

yx=[0;J'\b';0]+[ya;w;yb];  %matrix left division mldivide to obtain soln to linear system



end
%Test example on page 695 in chapter 11.4
%fun.f = @(x,y,dy) -(y./8).*dy+x.^3./4+4;
%fun.dfdy = @(x,y,dy) -(y./8);
%fun.dfy = @(x,y,dy) -(dy./8);
%fun.r = @(x,y,dy) x.^3./4+4;
%nlindiff( 1, 3, 17, 43/3, .1, fun, 10, 10^-5 )





%Now run the nonlinear BVP using nlinshoot and nlindiff
%fun.f = @(x,y,dy) 2*y^3-6*y-2*x^3;
%fun.dfdy = @(x,y,dy) 0;  %partial of f with respect to y'
%fun.dfy = @(x,y,dy) 6*y^2-6;  %partial of f with respect to y

%Exercise 11.2.4b
%[ yx2, error2 ] = bvpsolver( [1,2], [2,5/2], fun, .1, 'nlindiff' )