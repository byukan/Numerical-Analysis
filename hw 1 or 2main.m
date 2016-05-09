%main script solves IVPs specified in Solutions

%plug in start and end position, 
m = odesolver(0,20,[2,0],0.01,0.25,10^-4,@fvdp, 2,'rkf');
size = length(m);

for i = 1:size


ti(i) = m{i}(1);

x(i) = m{i}(2);

y(i) = m{i}(3);    

end


figure(1);

plot(ti,x,'-o')

figure(2);

plot(ti,y,'-o')

figure(3);

plot(x,y,'-o')

