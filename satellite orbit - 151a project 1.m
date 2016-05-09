function [root, a, b] = bis(a,b,tolb,func,pars)
%routine to find a root of a function using bisection method
a=a+pars;  %pars stores all constants in Kepler's equation
b=b+pars;
%func=@(x)(x-e*sin(x)-Mi);  %Kepler’s equation to find root E(t)

while abs(b-a)>tolb
p=a+(b-a)/2;  %compute midpoint p
	if (func(a)==0 || func(b)==0 || func(p)==0 || (b-a)/2<tolb)
    	root=p;
	end
    	if sign(func(a))*sign(func(p))>0  %same sign indicates root not in interval
                	a=p;  %consider second interval
    	else
            	b=p;  %consider first interval
    	end
end
root=p;
end

function [root, error, found, iter] = nwt(seed,toln,nmax,func,dfunc,pars)
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
	error=p-p0;
end
root=p;
found=false;  %the method failed
iter=i;  %iter is the number of iterations actually used
end

function [root, iter, itern] = bisnwt(a,b,tolb,toln,nmax,func,dfunc,pars)
%root finding algorithm using bisection method then using Newton's method
%to accelerate convergence
error=abs(b-a);  %absolute value of b-a is pi
while error>toln
	seed=bis(a,b,tolb,func,pars);  %first approximation to use in Newton
	[Ei,error]=nwt(seed,toln,nmax,func,dfunc,pars);
	if error>toln
    	tolb=tolb/2;
	else
    	break;
	end
end
root=Ei;
end

function [orbit] = tanom(T,e,n,varargin)
%application that solves the Kepler equation and gives the position of the
%satellite in its orbital plane
orbit=zeros(n+1,6);  %create matrix orbit to store values
for i=0:n
	ti=i*T/n;  %calculate time as a function of the orbit's period
	orbit(i+1,1)=ti;  %store value in matrix orbit
	Mi=2*pi*(ti/T);  %average anomaly M
	tolb= .1;  %error for bisection method
	a=0;
	b=pi;
	pars=0;
	if Mi>pi;
    	pars=pi;
	end
	toln=1e-16;  %suppose we want the error in the coordinates of the orbit
	%to be smaller than 10^-6
	nmax=100;  %test to validate code
	func=@(x)(x-e*sin(x)-Mi);  %Kepler’s equation to find root E(t)
	dfunc=@(x)1-e*cos(x);  %derivative of f
	Ei=bisnwt(a,b,tolb,toln,nmax,func,dfunc,pars);
	orbit(i+1,2)=Ei;  %use Ei to get true anomaly v
	vi=acos((e-cos(Ei))/(e*cos(Ei)-1));  %use v to get x,y position
	if Mi>pi;
    	vi=abs(vi-2*pi);
	end
	orbit(i+1,3)=vi;
	u=3.986012e5;  %constant mu to get a and r
	a=nthroot(u*((T*3600/(2*pi))^2), 3);
	ri=a*(1-e^2)/(1+e*cos(vi));
	orbit(i+1, 4)=ri;
	xi=ri*cos(vi);
	yi=ri*sin(vi);
	orbit(i+1 ,5)=xi;
	orbit(i+1 ,6)=yi;
end
end


