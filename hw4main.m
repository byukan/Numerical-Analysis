%First run the linear BVP using linshoot and lindiff
fun.f = @(x,y,dy) -(y./8).*dy+x.^3./4+4;
fun.dfdy = @(x)0;  %reduced to p(x) in the linear case
fun.dfy = @(x)-4;  %reduced to q(x) in the linear case
fun.r = @(x) cos(x);  %r(x) in the linear case

%Exercise 11.1.4b
[ yx1, error1, actual1 ] = bvpsolver( [0,pi/4], [0,0], fun, pi/20, 'linshoot' )
%Exercise 11.3.4b
[ yx3, error3, actual3 ] = bvpsolver( [0,pi/4], [0,0], fun, pi/20, 'lindiff' )




%Now run the nonlinear BVP using nlinshoot and nlindiff
fun.f = @(x,y,dy) 2*y^3-6*y-2*x^3;
fun.dfdy = @(x,y,dy) 0;  %partial of f with respect to y'
fun.dfy = @(x,y,dy) 6*y^2-6;  %partial of f with respect to y

%Exercise 11.2.4b
[ yx2, error2, actual2 ] = bvpsolver( [1,2], [2,5/2], fun, .1, 'nlinshoot' )
%Exercise 11.4.4b
[ yx4, error4,actual4 ] = bvpsolver( [1,2], [2,5/2], fun, .1, 'nlindiff' )



figure(1);
plot(yx1,'-ob')  %Exercise 11.1.4b
hold on
plot(actual1,'-ok');  %actual1 =actual3
title('Linear BVP Exercise 1');



figure(3);
plot(yx3, '-or')   %Exercise 11.3.4b
hold on 
plot(actual3,'-ok');  %actual1 =actual3
title('Linear BVP Exercise 3');



figure(2);
plot(yx2,'-ob')  %Exercise 11.2.4b
hold on
plot(actual2,'-ok');  %actual2=actual4
title('Non-Linear BVP Exercises 2');

figure(4);
plot(yx4, '-or')  %Exercise 11.4.4b
hold on
plot(actual4,'-ok');  %actual2=actual4
title('Non-Linear BVP Exercise 4');
