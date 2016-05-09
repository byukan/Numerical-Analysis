% Math 151A, Project 2 sample code

function S = spline3(t,ft,varargin)
% Computes the natural or cubic spline through data X,Y
% Algorithm modified slightly from Burden and Faires, 
% 'Numerical Analysis', Alg. 3.4 and 3.5
%
% Uses a directly implemented tridiagonal solve from the alg. noted above
%
% varargin{1} contains the spline type ('natural' or 'clamped')
% while varargin{2:3} contains the value of f' at the endpoints.
% Alternately, it can be passed in varargin{2} as a pair [f'(0) f'(x_n)].

%The spline is of the form
% S_j = f_j + b_j(x-x_j) + c_j(x-x_j)^2 + d_j(x-x_j)^3
%
% The first column of S is the function values f_j
% Columns 2-4 are b_j, c_j, and d_j respectively

n_vararg = length(varargin);
S(:,1) = ft;
S(:,2:4) = 0; 


%Process variable arguments to find spline type
s_type = 'natural'; %default
if(n_vararg>=1)
    s_type = varargin{1}; %override with input
    
    if(strcmp(s_type,'clamped'))
        if(n_vararg==3) % derivs. as two inputs
            [fpL, fpR] = deal(varargin{2:3});
        elseif(n_vararg==2) %derivs as a single array of two elements
            derivs = varargin{2};
            fpL = derivs(1);
            fpR = derivs(2);
        elseif(n_vararg==1 && length(t)>=5)
            %No derivative information given;
            %use a 5-point estimate to obtain boundary conditions
            fpL = deriv_fwd(t,ft);
            fpR = deriv_back(t,ft);            
        else 
            %If format is wrong, use natural spline
            fprintf('Wrong argument format. Defaulting to natural spline.\n');
            s_type = 'natural';
        end
    elseif(~strcmp(s_type,'natural'))
        %Only natural and clamped splines are supported
        fprintf('Unsupported type. Defaulting to natural spline.\n');
    end
end

n = length(t) - 1; %points are x_0 ... x_n
H = diff(t); %spacings

%Tridiagonal solve for Ac = b
%Find LU and solve Lz = b, with U unit upper-triangular
%Here mu is the upper diagonal of U
%and the vector L is the diagonal of the matrix L
%(the lower diagonal is just the lower diagonal of A)

L = zeros(1,n+1); 
Z = zeros(1,n+1);
mu = zeros(1,n+1);

%Do the first step explicitly (because of special cases for matrix entries)
if(strcmp(s_type,'clamped'))
    alphaL = 3*(S(2,1) - S(1,1))/H(1) - 3*fpL;
    alphaL = alphaL/(H(1)+H(2));
    L(1) = 2*H(1)/(H(1)+H(2));
    mu(1) = 0.5;
    Z(1) = alphaL/L(1);
else
    %L(1) = 1; 
    mu(1) = 0;
    Z(1) = 0;
end

%The remaining steps for finding L,U and Lz = b:
for i=1:n-1
    hr = H(i)/(H(i)+H(i+1));
    L(i+1) = 2 - hr*mu(i);
    mu(i+1) = (1-hr)/L(i+1); 
    alphai = 3*(ft(i+2)-ft(i+1))/H(i+1) - 3*(ft(i+1)-ft(i))/H(i);
    alphai = alphai/(H(i)+H(i+1));
    Z(i+1) = (alphai - hr*Z(i))/L(i+1); 
end

% Now solve Uc = z for c
% And also find the coefficients b,d

%Do the first step explicitly...
if(strcmp(s_type,'clamped'))
    alphaR = (3*fpR - 3*(S(n+1,1) - S(n,1))/H(n))/H(n);
    L(n+1) = 2-mu(n);
    Z(n+1) = (alphaR - Z(n))/L(n+1);
    S(n+1,3) = Z(n+1);
else %natural
    %L(n+1) = 1;
    Z(n+1) = 0;
    S(n+1,3) = 0;  
end

%...and then the remaining steps
for j=n-1:-1:0
    S(j+1,3) = Z(j+1) - mu(j+1)*S(j+2,3);
    S(j+1,2) = (S(j+2,1) - S(j+1,1))/H(j+1) - H(j+1)*(S(j+2,3) + 2*S(j+1,3))/3;
    S(j+1,4) = (S(j+2,3) - S(j+1,3))/(3*H(j+1));
end

S = S(1:n,:); %remove the last row

    %derivative estimates for clamped spline
    function dy = deriv_fwd(pts,vals)
        C = [-25/12, 4, -3, 4/3, -1/4];
        %assume equally spaced
        h = pts(2) - pts(1);
        dy = dot(C,vals(1:5))/h;
    end

    function dy = deriv_back(pts,vals)
        C = [1/4, -4/3, 3, -4, 25/12];
        %assume equally spaced
        h = pts(2) - pts(1);
        dy = dot(C,vals(end-4:end))/h;     
    end

end