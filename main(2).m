pc = odesolver(0,1,1,0.1,@ivp1,'euler');
size = length(pc);

for i = 1:size
t(i) = pc(i),(1);
w(i) = pc(i),(2); 
y(i) = sqrt(4-3*exp((-1)*t(i)^2));
d(i) = abs(y(i)-w(i));
end

figure(1);
plot(t,y,'-og',t,w,'-r')
%y is the actual soln in green circles, and w is the approximation in blue
%they overlap
title('Eulers method for ivp1(5.6.2d)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (2);  %plot of the global error
title('Global error between the approximated solution and the actual solution');
plot(d, '-b')



pc = odesolver(0,1,1,0.1,@ivp1,'am3');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = sqrt(4-3*exp((-1)*t(i)^2));
d(i) = abs(y(i)-w(i));
end

figure(3);
plot(t,y,'-og',t,w,'-r')
%y is the actual soln in green circles, and w is the approximation in blue
%they overlap
title('Adams-Moulton 3-step implicit method for ivp1(5.6.2d)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (4);  %plot of the global error
title('Global error between the approximated solution and the actual solution');
plot(d, '-b')




pc = odesolver(0,1,1,0.1,@ivp1,'adams4');
size = length(pc);

for i = 1:size
t(i) = pc{i}(1);
w(i) = pc{i}(2); 
y(i) = sqrt(4-3*exp((-1)*t(i)^2));
d(i) = abs(y(i)-w(i));
end

figure(3);
plot(t,y,'-og',t,w,'-r')
%y is the actual soln in green circles, and w is the approximation in blue
%they overlap
title('Adams-Moulton 3-step implicit method for ivp1(5.6.2d)');
xlabel('t');
ylabel('y(t) vs. w(t)');

figure (4);  %plot of the global error
title('Global error between the approximated solution and the actual solution');
plot(d, '-b')
