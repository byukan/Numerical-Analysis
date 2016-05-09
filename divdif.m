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