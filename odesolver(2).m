function [ yt ] = odesolver( a, b, y0, h, f, option )
%a is the initial time, b is the final time, y0 is the initial condition
%h is the step size by default, matlab array indicies start from 1, be
%aware when defining the size of yt using h
%f is the function handle f(t,y) by default
%option can be strings 'euler', 'rk4', 'adams4', 'am3'

switch option
    case 'euler';  %Algorithm 5.1, explicit one-step
        t=a;
        w=y0;
        N=(b-a)/h;
        
        for i=1:N
             t(i) = a+i*h;
             yt(i) = euler(t(i),y0,h,f);
             y0=yt(i);
         end
        
    case 'rk4';  %Algorithm 5.2, explicit one-step
        t=a;
        w=y0;
        yt = rk4( t, w, h, f );
        
    case 'adams4';  %Algorithm 5.4, explicit & implicit multistep
%updates wi using the Adams fourth order predictor-corrector method
        t=a;
        w=y0;
        a1=[55, -59, 37, -9]./24;
%a1 is the coefficient vector associated to wi, wi-1, wi-2, wi-3 in the prediction step
        a2=[9, 19, -5, 1]./24;
%a2 is the coefficient vector associated to wi+1, wi, wi-1, wi-2 in the correction step
        yt = adams4( t, w, h, f, a1, a2 );
        
        
    n=(b-a)/h;  %do iterate steps in odesolver since adams4.m does not have b
    for i=4:n
    t(i+1)=t(i)+h;
    yt(i+1) = yt(i)+h*(a1(1)*f(t(i),yt(i))+a1(2)*f(t(i-1),yt(i-1))+a1(3)*f(t(i-2),yt(i-2))+a1(4)*f(t(i-3),yt(i-3)));
    % Adams-Bashforth, predict wi
	yt(i+1) = yt(i)+h*(a2(1)*f(t(i+1),yt(i+1))+a2(2)*f(t(i),yt(i))+a2(3)*f(t(i-1),yt(i-1))+a2(4)*f(t(i-2),yt(i-2)));
    % Adams-Moulton, correct wi
    end
        
            
    case 'am3'; %implicit multistep
        t=a;
        w=y0;
        a2=[9, 19, -5, 1]./24;  %a2 is the coefficient vector associated to wi+1, wi, wi-1, wi-2
        n=(b-a)/h;
        
        yt{1,1}=[t,w];  %declare the solution matrix yt
        
        for i=1:2
            w = rk4(t,w,h,f);
            t = a+i*h;
            yt{i+1,1}=[t,w];  %update next row entry in yt
        end
        
        for j = 3:n

        t = [yt{j,1}(1),yt{j-1,1}(1),yt{j-2,1}(1)];  %t stores the 3 preceding times
        w = [yt{j,1}(2),yt{j-1,1}(2),yt{j-2,1}(2)];  %w stores the 3 preceding times
        
        yt{j+1,1}=am3(t, w, h, f, a2);  %update next time in yt
        end
    
end
        
end