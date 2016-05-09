function [ p ] = secant( p0, p1, f, tol, itermax )
%finds the root of the equation F(wi+1)=0 in the implicit multistep methods
%by using the secant method, see Algorithm 2.4
%p0 and p1 are the initial approximations
%TOL is the error tolerance
%itermax is the maximum number of iterations

i=2;
q0=f(p0);
q1=f(p1);

while i<=itermax
    p=p1-q1*(p1-p0)/(q1-q0);  %compute pi
    if abs(p-p1)<tol
        break;  %the procedure was successful
    end
    i=i+1;
    p0=p1;  %update p0, q0, p1, q1
    q0=q1;
    p1=p;
    q1=f(p);
    
end

end