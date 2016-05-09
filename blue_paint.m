function [blue_needed] = blue_paint(a,b,tolb,func,pars)
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