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
    vi=arccos((e-cos(Ei))/(e*cos(Ei)-1));  %use v to get x,y position
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