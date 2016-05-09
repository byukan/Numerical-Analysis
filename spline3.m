function S=spline3(t,ft,varargin) %computes the cubic splines corresponding
%to the interpolation points {(ti,f(ti))} and stores all coefficients in S
%to use CLAMPED, use format:  spline3(t,ft,fa,fb)

%construct NATURAL cubic spline interpolant S for the function f, defined at
%numbers xo,...,xn, satisfying S"(x0)=S"(xn)=0

if strcmp(varargin{1},'clamped')== 1;  %specify clamed in varargin
%construct CLAMPED cubic spline interpolant S for the function f defined at
%the numbers x.,...,xn, satisfying S'(x0)=f'(x0), and S'(xn)=f'(xn)
fa=varargin{2};  %contains f'(t0)
fb=varargin{3};  %contains f'(tn)
h=diff(t);
n=length(t)-1;
A=sparse(2:n+1,1:n,h,n+1,n+1) + ...
  sparse(1:n,2:n+1,h,n+1,n+1) + ...
  sparse(2:n,2:n,2*(h(1:n-1)+h(2:n)),n+1,n+1);
A(1,1)=2*h(1); A(n+1,n+1)=2*h(n);
b=[3./h(1)*(a(2)-a(1))-3*fa,3./h(2:n).*(a(3:n+1)-a(2:n))-3./h(1:n-1).*(a(2:n)-a(1:n-1)),3*fb-3/h(n)*(a(n+1)-a(n))]';
c=(A\b)';
b=(a(2:n+1)-a(1:n))./h-h./3.*(2*c(1:n)+c(2:n+1));
d=(c(2:n+1)-c(1:n))./(3*h);
c=c(1:n);
S=[ft',b',c',d'];  %store the vectors as columns in the matrix by using transpose

%by default when varargin is empty, we use the natural boundary condition
%without any declarations
else  %specity natural in varargin
h=diff(t);  %diff is the approximate derivative, h=x_i+1-x_i, which makes
%sense because it is the slope between two points.
n=length(t)-1;  %n is the length of the vector t
%construct the matrix A, which is a tri-diagonal matrix, using the sparse
%function from Matlab for efficiency.
A=sparse(2:n,1:n-1,h(1:n-1),n+1,n+1) + ...
  sparse(2:n,3:n+1,h(2:n),n+1,n+1) + ...
  sparse(2:n,2:n,2*(h(1:n-1)+h(2:n)),n+1,n+1);
A(1,1)=1;  A(n+1,n+1)=1;  %set the first and last diagonal entries to 1
%construct the splines
b=[0,3./h(2:n).*(a(3:n+1)-a(2:n))-3./h(1:n-1).*(a(2:n)-a(1:n-1)),0]';
c=(A\b)';
b=(a(2:n+1)-a(1:n))./h-h./3.*(2*c(1:n)+c(2:n+1));
d=(c(2:n+1)-c(1:n))./(3*h);
c=c(1:n);
%output S, which stores all coefficients, use a cell since they are
%different lengths
S=[ft',b',c',d'];  %store the vectors as columns in the matrix by using transpose
end