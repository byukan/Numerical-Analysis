function [ lambda, v, iter ] = invpower( A, x0, q, tol, maxiter)
%This function implements the inverse power method with l-2-norm scaling (refer to lecture notes)
%If the matrix A has eigenvalues lambda 1 to lambda n, and q is not an eigenvalue of A, then the matrix
%(A-q*I_n)^-1 has eigenvalues (lambda_i-q)^-1 with the same eigenvectors as A.  
%In other words, if mu is the dominant eigenvalue of (A-q*I_n)^-1, then mu^-1+q is the eigenvalue of A closest to q. 
%To obtain the dominant eigenvalue of A, you should choose the parameter q to be as close to the dominant
%eigenvalue of A as possible, that is not an eigenvalue of A.
%Notice that the power1 and power2 provide the approximate dominant eigenvalue of A.
%Inputs:
%A is the input matrix A.  In the PageRank algorithm, the modified adjacency matrix M should be plugged into A.
%x0, nonzero with unit length, is the initial guess for the dominant eigenvector
%q is a parameter such that lambda is the nearest eigenvalue to q
%tol is the tolerance
%maxiter is the maximal number of iterations allowed
%Outputs:
%lambda is the estimated dominant eigenvalue
%v is the estimated dominant eigenvector
%iter is the number of iterations that the algorithm actually used

k=1;
x=x0;
n=length(A);

while k<= maxiter
    
    x1=x;
    
    w=(A-q.*eye(n))\x;
    x=w/norm(w);
    lambda=x'*A*x;  %Rayleigh quotient

    if norm(x-x1)<tol
        break;
    end

    iter=k;

    k=k+1;
    
end

v=x;
end