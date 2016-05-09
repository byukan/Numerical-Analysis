function [root, found, iter] = nwt(seed,toln,nmax,func,dfunc,pars)
%routine to find root of a function using Newton's method
%found is a logical indicating whether root is found
i=1;
p0=seed;  %initial approximation
while i<=nmax
    p=p0-(func(p0)/dfunc(p0));
    if abs(p-p0)<toln
        found=true;
        root=p;  %the procedure was successful
        i=nmax;
    end
    i=i+1;
    p0=p;  %update p0
end
root=p;
found=false;  %the method failed
iter=i;  %iter is the number of iterations actually used
end