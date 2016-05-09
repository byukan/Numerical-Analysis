function [ wc ] = am3( t, w, h, f, a2)
%updates wi using implicit 3-step Adams-Moulton method
%t is a vector storing ti, ti-1, ti-2
%w is a vector storing wi, wi-1, wi-2
%a is the coefficient vector associated to wi+1, wi, wi-1, wi-2
%varargin is an optional input which can store the derivative function
%F'(w) if Newton's method is used
%since implicit method, apply root-finding method (secant) to find wi+1

tc=t(1)+h;  %tc is t(i+1)

g=@(x) -x+w(1)+h*(a2(1)*f(tc,x)+a2(2)*f(t(1),w(1))+a2(3)*f(t(2),w(2))+a2(4)*f(t(3),w(3)));
%set anonymous function, x is yt(i+1)
yt_iplus1=secant( w(1), w(2), g, 10^-4, 100);  %pass into secant.m to find root

wc=[tc,yt_iplus1];

end