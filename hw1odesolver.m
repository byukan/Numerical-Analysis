function [yt] = hw1odesolver(a, b, y0, N, f, argf, option)

% initial time a, final time b, initial condition y0
% number of steps N (b-a=N*h), function that defines the IVP f
% argf is a vector stacking all parameters required by f, set argf=[] if no
% parameters are involved in f
% option is a string specifying the method being used, which can be 'euler',
%'midpoint', or 'meuler', use command switch..case
% yt is a vector containigng all approximations of the solution at given
% times ti

h=(b-a)/N;  % define the step size
y_initial=y0;
switch option
    case 'euler';
for i=1:N
t(i) = a+i*h;
yt(i) = euler(t(i),y_initial,h,f,argf);
y_initial=yt(i);
end
    case 'midpoint'
for i=1:N
t(i) = a+i*h;
yt(i) = midpoint(t(i),y_initial,h,f,argf);

y_initial=yt(i);
end
    case 'meuler'
for i=1:N
t(i) = a+i*h;
yt(i) = meuler(t(i),y_initial,h,f,argf);

y_initial=yt(i);
end
end
end







