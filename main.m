%solves the specified IVPS in solutions, exercise 5.6.2d, and exercise
%5.11.2b with h=.1, and h=2

pc = odesolver(0,1,1,0.1,@ivp1,'euler');
size = length(pc);

for i = 1:size
t(i) = pc(i),(1);
w(i) = pc(i),(2); 
y(i) = sqrt(4-3*exp((-1)*t(i)^2));
d(i) = abs(y(i)-w(i));
end

figure(1);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Eulers method for ivp1(5.6.2d)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (2);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');



pc = odesolver(0,1,1,0.1,@ivp1,'rk4');
size = length(pc);

for i = 1:size
t(i) = pc(i),(1);
w(i) = pc(i),(2); 
y(i) = sqrt(4-3*exp((-1)*t(i)^2));
d(i) = abs(y(i)-w(i));
end

figure(3);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Runge-Kutta order 4 method for ivp1(5.6.2d)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (4);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');




pc = odesolver(0,1,1,0.1,@ivp1,'am3');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = sqrt(4-3*exp((-1)*t(i)^2));
d(i) = abs(y(i)-w(i));
end

figure(5);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Adams-Moulton 3-step implicit method for ivp1(5.6.2d)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (6);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');





pc = odesolver(0,1,1,0.1,@ivp1,'adams4');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = sqrt(4-3*exp((-1)*t(i)^2));
d(i) = abs(y(i)-w(i));
end

figure(7);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Adams fourth order predictor-corrector method for ivp1(5.6.2d)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (8);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');






%======Now for ivp2, exercise 5.11.2b with h=.1==========

pc = odesolver(0,1,exp(1),0.1,@ivp2,'euler');
size = length(pc);

for i = 1:size
t(i) = pc(i),(1);
w(i) = pc(i),(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(9);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Eulers method for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (10);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');




pc = odesolver(0,1,exp(1),0.1,@ivp2,'rk4');
size = length(pc);

for i = 1:size
t(i) = pc(i),(1);
w(i) = pc(i),(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(11);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('RK4 for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (12);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');



pc = odesolver(0,1,exp(1),0.1,@ivp2,'am3');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(13);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('3 step Adams-Moulton Implicit for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (14);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');



pc = odesolver(0,1,exp(1),0.1,@ivp2,'adams4');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(15);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Adams 4th order predictor-corrector for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (16);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');





%======Now for ivp2, exercise 5.11.2b with h=.2==========

pc = odesolver(0,1,exp(1),0.2,@ivp2,'euler');
size = length(pc);

for i = 1:size
t(i) = pc(i),(1);
w(i) = pc(i),(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(17);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Eulers method for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (18);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');




pc = odesolver(0,1,exp(1),0.2,@ivp2,'rk4');
size = length(pc);

for i = 1:size
t(i) = pc(i),(1);
w(i) = pc(i),(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(19);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('RK4 for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (20);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');



pc = odesolver(0,1,exp(1),0.2,@ivp2,'am3');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(21);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('3 step Adams-Moulton Implicit for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (22);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');



pc = odesolver(0,1,exp(1),0.2,@ivp2,'adams4');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = exp(-10*t(i)+1)+t(i);
d(i) = abs(y(i)-w(i));
end

figure(23);
plot(t,y,'-ok',t,w,'-r')
%y is the actual soln in black circles, and w is the approximation in blue
%they overlap
title('Adams 4th order predictor-corrector for ivp2(5.11.2b)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (24);  %plot of the global error
plot(d, '-b')
title('Global error between the approximated solution and the actual solution');


