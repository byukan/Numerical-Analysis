%solves the nonlinear problems specified in "solutions"
fun.F = @(x)[15*x(1)+x(2)^2-4*x(3)-13 ; x(1)^2+10*x(2)-x(3)-11 ; x(2)^3-25*x(3)+22];
fun.J = @(x) [15 2*x(2) -4 ; 2*x(1) 10 -1; 0 2*x(2)^2 -25];
paras.tol=10^-6;
paras.maxiter=10;

%Exercise 10.2.6c
paras.option='newton';
fsolver( [0;0;0], fun, paras )
%Exercise 10.3.5c
paras.option='broyden';
fsolver( [0;0;0], fun, paras )
%Exercise 10.4.2a
paras.option='steep';
fsolver( [0;0;0], fun, paras )
%Exercise 10.5.4c
paras.option='homotopy';
fsolver( [0;0;0], fun, paras )



figure(1);
plot(newton_x,'-ob')  %Exercise 10.2.6c
title('Newton');


%odesolver( 0, 1, [0;0;0], 1/10, fun, [], 'rk4')



%[yt] = hw1odesolver(0, 1, [0;0;0], 10, fun, [], 'midpoint')