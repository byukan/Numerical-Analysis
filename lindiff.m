function [ yx ] = lindiff( a, b, ya, yb, h, p, q, r )
%applies the finite-difference method to solve the linear BVP (2)
%In order to solve the tridiagonal system of the form Aw=b, I first
%construct the coefficient matrix A by using the command
%spdiags([d1,d2,d3], -1:1,N,N) which creates a NxN sparse band matrix.
%Then I use the matrix left division mldivide, ie A/b, to obtain the
%solution to the linear system
%a and b are endpoints of BVP
%ya is boundary condition alpha
%yb is boundary condition beta
%h is the step size
%p,q,r are the respective function handles of p(x), q(x) and r(x) of BVP
%output yx is a vector storing the approximations of the solution at xi

n=(b-a)/h;  %difference methods step size, use same formula for this algorithm


for i=1:n-1  %sparse, spdiags will automatically cut off d1(1) and d3(n), that's how it
%works, so make all the diagonals the same size.
    x(i)=a+i*h;  %x is the vector of time steps used to calculate entries of d1, d2, and d3
    d1(i)=-1-(h/2).*p(x(i));  %d1 is the lower diagonal of tridiagonal matrix
    d3(i)=-1+(h/2).*p(x(i));  %d3 is the upper diagonal of tridiagonal matrix
    d2(i)=2+h^2.*q(x(i));  %d2 is the main diagonal of tridiagonal matrix
end


b(1)=-h^2.*r(x(1))+(1+(h/2).*p(x(1))).*ya;  %%construct vector b, starting with first entry
b(n-1)=-h^2.*r(x(n-1))+(1-(h/2).*p(x(n-1))).*yb;  %last entry
for i=2:n-2  %fill in inner entries of b
    b(i)=-h^2.*r(x(i));
end

A=spdiags([d3',d2',d1'], -1:1, n-1, n-1);  %A is the coefficient matrix
%to make spdiags work, I had to take the transpose so the size works with
%this function, also this function is written to use n-1 as size
yx=[ya;A'\b';yb];  %matrix left division mldivide to obtain soln to linear system

end
%Test example on page 688 in chapter 11.3
%linshoot(1,2,1,2, .1, @(x)(-2/x), @(x)(2/x^2), @(x) sin(log(x))/x^2)