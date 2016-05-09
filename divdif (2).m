% Math 151A, Project 2 sample code

function [P, Q] = divdif(X,Y)
% Taken from Algorithm 3.2 of Burden and Faires
%
% Computes the Newton divided difference table for data (X,Y)
% and returns the coefficients of the interpolating polynomial P_n
% in the form c_(n-1), ... , c_0
%
% Inputs:
%   X - the x-values
%   Y - the function values
%
% Outputs:
%   P - the polynomial coefficients (highest degree to lowest)
%   Q - the divided difference table, in lower triangular form 

n = length(X); % polynomial degree + 1

if(nargout > 1)
    %Compute the full matrix of coefficients
    Q(:,1) = Y;
    Q(:,2:n) = 0;

    for i=1:n-1
        for j=i+1:n
            Q(j,i+1) = (Q(j,i) - Q(j-1,i))/(X(j) - X(j-i));
        end
    end
    P = diag(Q);
    
else
    %Just compute the polynomial coefficients
    for i=1:n-1
       for j=n:-1:i+1
           Y(j) = (Y(j) - Y(j-1))/(X(j) - X(j-i));
       end
    end
    
    P = Y; 
end