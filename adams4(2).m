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

a=t;
yt=w;

    for i=1:3  %compute starting values using Runge-Kutta method
     k1=h*f(t(i),yt(i,:));
     k2=h*f(t(i)+h/2,yt(i,:)+k1/2);
     k3=h*f(t(i)+h/2,yt(i,:)+k2/2);
     k4=h*f(t(i)+h,yt(i,:)+k3);
     yt(i+1,:)=yt(i,:)+(k1+2*k2+2*k3+k4)/6;
     t(i+1)=t(i)+h;
    end
    
wc=yt;
    
end

