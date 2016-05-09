% Math 151A, Project 2 sample code
%
% main function for generating requested project results
% See the problem statement for details
%
% NOTE: spline3_spdiags uses MATLAB's sparse matrix type
% to do the tridiagonal solve, while spline3 directly implements
% an LU factorization.
%
% There are a few more plots than required; they are here to illustrate
% interesting features or to show other ways to plot the data.
%
% I cannot guarantee this code is bug free or conforms to proper style.

%Load the given data
data = load('data.mat');
ip = data.ip;
tp = data.tp;
tu = linspace(0,3.1,311)';

%Construct the paths
[path_lagrange, L] = findpath(ip,tp,'lagrange');
[path_spline, S] = findpath(ip,tp,'spline');
[path_lagrange_u, L_u] = findpath(ip,tu,'lagrange');
[path_spline_u, S_u] = findpath(ip,tu,'spline');

%save the data (change names to agree with problem statment)
path1a = path_lagrange;
path1b = path_spline;
path2a = path_lagrange_u;
path2b = path_spline_u;

save('results.mat','path1a','path1b','path2a','path2b');

%Plot the results
figure(1)
title('random sample times')
subplot(1,2,1)
plot(path_lagrange(:,2),path_lagrange(:,3),'-k', ...
    ip(:,2),ip(:,3), '--.r');
axis([-1 1 -1 1.3]);
xlabel('x');
ylabel('y');
%title('Lagrange interpolant');
legend('Lagrange','data');

subplot(1,2,2)
plot(path_spline(:,2),path_spline(:,3),'-k', ...
    ip(:,2),ip(:,3), '--.r');
axis([-1 1 -1 1.3]);
xlabel('x');
ylabel('y');
%title('cubic spline');
legend('spline','data');

set(gcf,'PaperPositionMode','auto');
set(gcf,'Units','inches');
set(gcf,'OuterPosition',[2 2 8 4]);

%An alternate plot
figure(2)
title('uniform sample times')
subplot(1,2,1)
plot(tu,path_lagrange_u(:,2),'-k',ip(:,1),ip(:,2),'.r',...
    tu,path_lagrange_u(:,3),'--k',ip(:,1),ip(:,3), '*b');
axis([-0.5 3.1 -1.2 1]);
xlabel('t');
ylabel('x,y');
%title('Lagrange interpolant');
legend('Lagrange, x','Lagrange, y','data, x','data, y');

subplot(1,2,2)
plot(tu,path_spline_u(:,2),'-k',ip(:,1),ip(:,2),'.r',...
    tu,path_spline_u(:,3),'--k',ip(:,1),ip(:,3), '*b');
axis([-0.5 3.1 -1.2 1]);
xlabel('t');
ylabel('x,y');
%title('Lagrange interpolant');
legend('Spline, x','Spline, y','data, x','data, y');

%Needed to ensure that 'print' produces the 
%correct dimensions for the plot
set(gcf,'PaperPositionMode','auto');

%Set the plot window size
set(gcf,'Units','inches');
set(gcf,'OuterPosition',[2 2 8 4]);