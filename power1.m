function [ lambda, v, iter ] = power1( A, x0, tol, maxiter)
%This function implements the power method with l-infinity-norm scaling (refer to algorithm 9.1)
%Note that in algorithm 9.1, the Rayleigh quotient is replaced by the pth component 
%of the vector Ax^(k-1) where p is the smallest integer in [1,n] such that abs(x_p) = norm(x, Inf).
%Inputs:
%A is the input matrix A.  In the PageRank algorithm, the modified adjacency matrix M should be plugged into A.
%x0 is the initial guess for the dominant eigenvector
%tol is the tolerance
%maxiter is the maximal number of iterations allowed
%Outputs:
%lambda is the estimated dominant eigenvalue
%v is the estimated dominant eigenvector
%iter is the number of iterations that the algorithm actually used
%For faster convergence, I use the method presented in lecture using
%Rayleigh quotient

k=1;
x=x0;  %initialize guess for dominant eigenvector

while k<=maxiter
    x=A*x;  %apply A
    lambda=(x'*A*x)/(x'*x);
    
    x1=x/norm(x,Inf);
    if norm(x-x1, Inf)<tol
        x=x1;
        break;
    end
    x = x1;  %update x

    iter=k;

    k=k+1;  %the drawback of this version is lim as k tends to Inf of lambda1^k=Inf if lambda1>1

end
v=x;
end