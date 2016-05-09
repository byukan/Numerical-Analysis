% Jeffrey Wong, 3/11/14
% Math 151A, Project 2 sample code

function S = spline3_spdiags(t,ft,varargin)
% Computes the natural or cubic spline through data X,Y
% Algorithm modified slightly from Burden and Faires, 
% 'Numerical Analysis', Alg. 3.4 and 3.5
%
% Uses spdiags to build the tridiagonal matrix and backslash to solve
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

%Space for the spline coefficients
S(:,1) = ft(1:n);

H = diff(t); %spacings; H(i+1) = x_(i+1) - x_i, i=0,..n-1 
mu(2:n) = H(2:n)./(H(1:n-1)+H(2:n)); %superdiagonal
lambda(1:n-1) = 1 - mu(2:n); %sub diagonal
dY = diff(ft); %dY(i+1) = y(i+1) - y(i), i = 0,..,n-1

%An incredibly obtuse way of finding the RHS
r(2:n) = 6*(dY(2:n)./H(2:n) - dY(1:n-1)./H(1:n-1))./(H(1:n-1)+H(2:n));

if(strcmp(s_type,'natural'))
   mu(1) = 0;
   lambda(n) = 0;
   r(1) = 0;
   r(n+1) = 0;
else
    mu(1) = 1;
    lambda(n) = 1;
    r(1) = 6/H(1)^2*(dY(1)) - 6/H(1)*fpL;
    r(n+1) = 6/H(1)*fpR - 6*dY(n)/H(n)^2;
end

%needed for spdiags
lambda(n+1) = 0; 
mu = [0 mu]; 

A = spdiags([lambda' 2*ones(n+1,1) mu'],[-1 0 1],n+1,n+1);
m = (A\r')'; %solve the system, save as row vector
         %Note that because A is sparse, it will use a sparse LU algorithm.
b = dY(1:n)./H(1:n) - H(1:n).*(m(2:n+1) + 2*m(1:n))/6;
c = m(1:n)/2;
d = 1/6*diff(m)./H(1:n);

S = [S b' c' d'];

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