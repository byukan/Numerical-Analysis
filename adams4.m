function [ wc ] = adams4( t, w, h, f, a1, a2 )
%part of Algorithm 5.4, explicit & implicit multistep
%updates wi using the Adams fourth order predictor-corrector method
%t is a vector which stores the four preceding times ti, ti-1, ti-2, ti-3
%to generate wi+1
%w is a vector storing wi, wi-1, wi-2, wi-3
%a1 is the coefficient vector associated to wi, wi-1, wi-2, wi-3 in the
%prediction step, ie [55, -59, 37, -9]./24, which can be defined in
%odesolver.m
%a2 is the coefficient vector associated to wi+1, wi, wi-1, wi-2 in the
%correction step, ie [9, 19, -5, 1]./24, which can be defined in
%odesolver.m
%wc is the generated wi+1, which is a scalar in first-order ODEs but a
%vector in systems of ODEs, use vector case for more general code

tc=t(1)+h;  %tc is t(i+1)
predictor = w(1) + h.*(a1(1).*f(t(1),w(1))+ a1(2).*f(t(2),w(2)) + a1(3).*f(t(3),w(3)) + a1(4).*f(t(4),w(4))); 
    % Adams-Bashforth, predict wi
corrector = w(1) + h.*(a2(1).*f(tc,predictor)+ a2(2).*f(t(1),w(1)) + a2(3).*f(t(2),w(2)) + a2(4).*f(t(3),w(3)));
    % Adams-Moulton, correct wi
    
wc=[tc, corrector];
    
end

