function [ lambda, v, iter ] = power2( A, x0, tol, maxiter)
%This function implements the power method with l-2-norm scaling (refer to algorithm 9.2 or lecture notes)
%Inputs:
%A is the input matrix A.  In the PageRank algorithm, the modified adjacency matrix M should be plugged into A.
%x0 is the initial guess for the dominant eigenvector
%tol is the tolerance
%maxiter is the maximal number of iterations allowed
%Outputs:
%lambda is the estimated dominant eigenvalue
%v is the estimated dominant eigenvector
%iter is the number of iterations that the algorithm actually used

k=1;
x=x0;  %initialize guess for dominant eigenvector

while k<=maxiter
    x=A*x;  %apply A
    lambda=(x'*A*x)/(x'*x);
    
    x1=x/norm(x,2);
    if norm(x-x1, Inf)<tol
        x=x1;
        break;
    end
    x = x1;  %update x

    iter=k;

    k=k+1;

end
v=x;

end