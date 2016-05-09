function [ lambda, v, iter ] = eigfinder( A, x0, paras )
%To find the dominant eigenvector and eigenvalue of a matrix, we define these MATLAB functions.  
%eigfinder.m finds the dominant eigenvector and teh associated eigenvalue of a given matrix A in R^nxn.
%Inputs:
%A is the input matrix A.  In the PageRank algorithm, the modified adjacency matrix M should be plugged into A.
%x0 is the initial guess for the dominant eigenvector
%paras is a structure array containing the following fields:
%- .tol is the tolerance
%- .maxiter is the maximal number of iterations allowed
%- .q is the parameter q used in the inverse power method and is empty in other methods
%- .option specifies the method to be used, eg., strings ‘power1’, ‘power2’, and ‘invpower’
%Outputs:
%lambda is the estimated dominant eigenvalue
%v is the estimated dominant eigenvector
%iter is the number of iterations that the algorithm actually used

switch paras.option
    case 'power1'
        
        [ lambda, v, iter ] = power1( A, x0, paras.tol, paras.maxiter);
        
    case 'power2'
        
        [ lambda, v, iter ] = power2( A, x0, paras.tol, paras.maxiter);
        
    case 'invpower'
        
        [ lambda, v, iter ] = invpower( A, x0, paras.q, paras.tol, paras.maxiter);

end

end