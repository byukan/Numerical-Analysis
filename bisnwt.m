function [root, iter, itern] = bisnwt(a,b,tolb,toln,nmax,func,dfunc,pars)
%root finding algorithm using bisection method then using Newton's method
%to accelerate convergence
error=pi;  %absolute value of b-a is pi
while error>toln
    seed=bis(a,b,tolb,func,pars);
    [Ei, found, itern]=nwt(seed,toln,nmax,func,dfunc,pars);
    if error>toln
        tolb=tolb/2;
    else
        break;
    end
end
root=Ei;
end