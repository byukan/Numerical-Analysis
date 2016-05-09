function [P, Q]=divdif(x,y)  %P stores all coefficients of Lagrange polynomial
%Q is a matrix with entries Fi,j=[xi-j,...,xj]
%implements algorithm 3.2 in the textbook
%Newton's Divided-Difference Formula
n=length(x)-1;  %n is length of x-1
F=zeros(n+1,n+1);  %F is square matrix of size n+1
F(:,1)=y(:);  %input y gets assigned to columns of F
for i=1:n
    for j=1:i
        F(i+1,j+1)=(F(i+1,j)-F(i,j))/(x(i+1)-x(i-+1));
    end
end
P=zeros(n+1); %put diagonal entries of F into vector P
    for i=1:n+1;
        P(i)=F(i,i);
    end
Q=F;
end


function pt=evalP(P,t_int,t_test)
%evaluates the polynomial with coefficients P associated with interpolation
%points t_int at the test points t_test and returns that value in pt.
%task is to evaluate the polynomial P(x)=(x-x0)*Q(x)+b0 at x0.
%work recursively backwards, starting at Fnn(x-x_n-1)
n=length(P)-1;  %gives n, the degree of P
pt=P(n+1)*(t_test-t_int(n));  %compute/starting at Fnn(x-x_n-1)
for j=n:2  %working backwards, using nested multiplication
    pt=(pt+P(j))*(t_test-t_int(j-1));  %compute bj for P
end
pt=pt+P(1);  %output pt, that last term is F(0,0)
end



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


function st=evalS(S,t_int,t_test)
%one method to evaluate the spline with coefficients stored in S at the
%point(s) t_test and returns that value in st.
n=length(t_test)-1;
st=0*t_int;
x=t_test;  xx=t_int;  %declare the test points and interpolation points
a=S(:,1);  %fist column of S gets assigned to a
b=S(:,2);  %second column of S gets assigned to b
c=S(:,3);  %third column of S gets assigned to c
d=S(:,4);  %fourth column of S gets assigned to d
for i=1:n
  ix=xx>=x(i) & xx<=x(i+1);
  st(ix)=a(i)+b(i)*(xx(ix)-x(i))+c(i)*(xx(ix)-x(i)).^2+d(i)*(xx(ix)-x(i)).^3;
end



function [path, coeff]=findpath(ip,tp,method,varargin)
%finds the path matrix of the form (1) and the coefficient matrix of the
%form (2) or (4).  Here ip is a matrix of the form (1) storing all
%interpolation points, tp is a column vector storing t0,...tm-1, and method
%specifies the interpolation method, eg. strings "lagrange" and "spline"

if strcmp(method,'lagrange')== 1;  %specify lagrange for interpolation method
    %construct & evaluate a single polynomial that interpolates all n nodes
    %first find coefficients of lagrange polynomial that gives x coordinates
    
     [Px,Qx]=divdif(ip(:,1),ip(:,2));  %use columns 1(time) and 2 x(t) of ip, which stores interpolation points
         x_coord=evalP(Px,ip(:,1),tp);  %Px gives x coordinate for each time t

     [Py,Qy]=divdif(ip(:,1),ip(:,3));  %use columns 1(time) and 3 y(t) of ip, which stores interpolation points
         y_coord=evalP(Py,ip(:,1),tp);  %Py gives y coordinate for each time t

            path=[ip(:,1),x_coord,y_coord];  %path matrix that gives the x and y coordinates at each time t
            coeff=[ip(:,1),Px,Py];  %save the polynomial coefficients
end





if strcmp(method,'spline')== 1;  %specify spline for interpolation method
   Sx=spline3(ip(:,1),ip(:,2),varargin); 
   %use spline3 to evaluate Sx using given interpolation points t, and x(t)
   Sy=spline3(ip(:,1),ip(:,3),varargin);  
   %use spline3 to evaluate Sy using given interpolation points t, and y(t)
   
   sx_coord=evalS(Sx,ip(:,1),tp);  %evaluate spline w/ coefficients stored in S
   sy_coord=evalS(Sy,ip(:,1),tp);  %at given points tp, and returns sx,sy coordinates
   
   path=[ip(:,1), sx_coord, sy_coord]; %path matrix that gives x,y coords at each time t
   coeff=[ip(:,1),Sx(:,1),Sx(:,2),Sx(:,3),Sx(:,4),Sy(:,1),Sy(:,2),Sy(:,3),Sy(:,4) ];
   %final matrix of the form (4), save all coefficients
   
end

