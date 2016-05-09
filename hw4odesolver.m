function [ yt ] = odesolver( a, b, y0, h, fun, argf, option )
%solves the system of first-order IVPs using one-step methods or multistep
%methods.  Notice that f: [a,b]XR^k--->R^k is a vector-valued function and
%yt returns a (N+1) X k matrix.  At least, you should include the
%Runge-kutta method of order four, and define a nested function rk4 in the
%odesolver.m.  Note that argf will be very useful when solving the z(x) -
%involved IVP in the nonlinear shooting method (z(x) =z(x,t)) ie
%z''(x)=fy(x,y,y')*z(x) + fy'(x,y,y')*z'(x), where y,y' are from the
%y(x)-involved IVP and will be updated for each tk.  You could use this
%routine to get a solution to a high-order IVP

%a is the initial time, b is the final time, y0 is the initial condition
%h is the step size
%f is the function fun, a structure array containing the fields .f, .dfdy,
%.dfy, and .r
%argf contains related parameters for fun, which is optional
%option can specify 'linshoot', 'nlinshoot', 'lindiff', and 'nlindiff'

switch option
    case 'rk4linshoot'  %linear shooting method
        n = (b-a)/h;  % step 1, n is number of subintervals
        u1=y0;
        u2=0;
        v1=0;
        v2=1;
        u = zeros(2,n);
        v = zeros(2,n);
        
        p=fun.dfdy;  %rename to use equations
        q=fun.dfy;
        r=fun.r;

    for i = 1 : n  % step 2, do steps 3 and 4, use Runge-Kutta method for systems

        x = a+(i-1)*h;  % step 3
        t = x+0.5*h;

        k11=h*u2;  % step 4
        k12=h*(p(x)*u2+q(x)*u1+r(x));
        k21=h*(u2+0.5*k12);
        k22=h*(p(t)*(u2+0.5*k12)+q(t)*(u1+0.5*k11)+r(t));
        k31=h*(u2+0.5*k22);
        k32=h*(p(t)*(u2+0.5*k22)+q(t)*(u1+0.5*k21)+r(t));
        t=x+h;
        k41=h*(u2+k32);
        k42=h*(p(t)*(u2+k32)+q(t)*(u1+k31)+r(t));
        u1=u1+(k11+2*(k21+k31)+k41)/6;
        u2=u2+(k12+2*(k22+k32)+k42)/6;
        k11=h*v2;
        k12=h*(p(x)*v2+q(x)*v1);
        t=x+0.5*h;
        k21=h*(v2+0.5*k12);
        k22=h*(p(t)*(v2+0.5*k12)+q(t)*(v1+0.5*k11));
        k31=h*(v2+0.5*k22);
        k32=h*(p(t)*(v2+0.5*k22)+q(t)*(v1+0.5*k21));
        t=x+h;
        k41=h*(v2+k32);
        k42=h*(p(t)*(v2+k32)+q(t)*(v1+k31));
        v1=v1+(k11+2*(k21+k31)+k41)/6;
        v2=v2+(k12+2*(k22+k32)+k42)/6;

        u(1,i) = u1;
        u(2,i) = u2;
        v(1,i) = v1;
        v(2,i) = v2;
    end
    
    
    yt{1}=u;  %make yt a cell to store vectors u and v to get passed back to linshoot.m
    yt{2}=v;
    
    
    
    
    case 'rk4nlinshoot'  %nonlinear shooting method
        
        w1=argf{1};
        w2=argf{2};
        
        u1 = 0;
        u2 = 1;
    
        n = (b-a)/h;  % step 1, n is number of subintervals

     for i = 1 : n   % step 4, the Runge-Kutta method for systems is used in steps 5 and 6
     
       x = a+(i-1)*h;  %  step 5
       T = x+0.5*h;

       k11 = h*w2(i);  % step 6
       k12 = h*fun.f(x,w1(i),w2(i));
       k21 = h*(w2(i)+0.5*k12);
       k22 = h*fun.f(T,w1(i)+0.5*k11,w2(i)+0.5*k12);
       k31 = h*(w2(i)+0.5*k22);
       k32 = h*fun.f(T,w1(i)+0.5*k21,w2(i)+0.5*k22);
       k41 = h*(w2(i)+k32);
       k42 = h*fun.f(x+h,w1(i)+k31,w2(i)+k32);
       w1(i+1) = w1(i)+(k11+2*(k21+k31)+k41)/6;
       w2(i+1) = w2(i)+(k12+2*(k22+k32)+k42)/6;
       k11 = h*u2;
       k12 = h*(fun.dfy(x,w1(i),w2(i))*u1+fun.dfdy(x,w1(i),w2(i))*u2);
       k21 = h*(u2+0.5*k12);
       k22 = h*(fun.dfy(T,w1(i),w2(i))*(u1+0.5*k11)+fun.dfdy(T,w1(i),w2(i))*(u2+0.5*k21));
       k31 = h*(u2+0.5*k22);
       k32 = h*(fun.dfy(T,w1(i),w2(i))*(u1+0.5*k21)+fun.dfdy(T,w1(i),w2(i))*(u2+0.5*k22));
       k41 = h*(u2+k32);
       k42 = h*(fun.dfy(x+h,w1(i),w2(i))*(u1+k31)+fun.dfdy(x+h,w1(i),w2(i))*(u2+k32));
       u1 = u1+(k11+2*(k21+k31)+k41)/6;
       u2 = u2+(k12+2*(k22+k32)+k42)/6;
       
     end
    
     
    yt{1}=w1;  %make yt a cell to store vectors u and v to get passed back to linshoot.m
    yt{2}=u1;        
        
        
        
        

end


end