% Jeffrey Wong, 3/11/14
% Math 151A, Project 2 sample code

%Script to demonstrate the Runge phenomenon
%Interpolates using Lagrange (equally spaced points)
%and Lagrange using the Chebyschev nodes.
%
%Also uses spline interpolation. 
%The plotting here is much more elaborate than necessary.

fR = @(x) 1./(1+25*x.^2);
fp = @(x) -((50*x)./(1+25*x.^2).^2);

Ns = [10 20 50];
N_tests = length(Ns);
X_sample = linspace(-1,1,1000)'; %Sample points to test the interpolant

figure(2)
clf
styles = {'-k','-r','-b'};
err_c = zeros(N_tests,1); 
err_spline = zeros(N_tests,1);
err_lagrange = zeros(N_tests,1);
leg_str = {}; %The plot legend, to be defined

for i=1:N_tests

%Define the nodes and function values at nodes
N = Ns(i);
nodes_equal = linspace(-1,1,N)';
nodes_ch = cos((2*(1:N)-1)*pi/(2*N))';

Y = fR(nodes_equal);
Yc = fR(nodes_ch);

%Build the interpolants
P_lagrange = divdif(nodes_equal,Y);
P_cheby = divdif(nodes_ch,Yc);
S_clamped = spline3(nodes_equal,Y,'clamped',fp(-1),fp(1));

%Evaluate at the test points for each interpolant
Y_spline = evalS(S_clamped, nodes_equal, X_sample);
Y_lagrange = evalP(P_lagrange,nodes_equal,X_sample);
Y_cheby = evalP(P_cheby,nodes_ch,X_sample);
Y_actual = fR(X_sample);

%Compute the error
err_lagrange(i) = max(abs(Y_lagrange - Y_actual));
err_c(i) = max(abs(Y_cheby - Y_actual));
err_spline(i) = max(abs(Y_spline - Y_actual));

%Add the appropriate interpolants to each subplot
figure(2)
subplot(1,3,1)
hold on %have plot add new lines, instead of overwriting
plot(X_sample,Y_lagrange,styles{mod(i,3)+1});

%Save the legend entry to be applied after the loop
leg_str{i} = ['N=',num2str(N)]; %#ok<SAGROW>

figure(2)
subplot(1,3,2)
hold on
plot(X_sample,Y_cheby,styles{mod(i,3)+1});

subplot(1,3,3)
hold on
plot(X_sample,Y_spline,styles{mod(i,3)+1});
end

%------Plotting the approximations---------
%Decorate each plot with labels, etc.
figure(2)
subplot(1,3,1)
plot(X_sample,fR(X_sample),'--k');
xlabel('x');
ylabel('P_n(x)');
legend([leg_str {'actual'}]); %add the legend
hold off

subplot(1,3,2)
plot(X_sample,fR(X_sample),'--k');
xlabel('x');
ylabel('P_n(x) (cheby.)');
hold off

subplot(1,3,3)
plot(X_sample,fR(X_sample),'--k');
xlabel('x');
ylabel('S(x)');
hold off

%This is needed for printing the plot via 'print',
%otherwise it won't use the correct dimensions
set(gcf,'PaperPositionMode','auto');

%Set the size of the plot window
set(gcf,'Units','inches');
set(gcf,'OuterPosition',[2 2 9 4]);

%-------Plotting the error----------
figure(3)
subplot(1,2,1)
semilogy(Ns,err_lagrange,'.-b');
xlabel('N');
ylabel('abs. error (L^\infty)');
axis([0 50 1e-1 1e5]);
legend('equally spaced');
%title('Absolute error for equally spaced Lagrange interpolation')

subplot(1,2,2)
semilogy(Ns,err_c,'.-k',Ns,err_spline,'.-r');
xlabel('N');
ylabel('abs. error (L^\infty)');
legend('chebyschev','spline');
axis([0 50 1e-4 1e0]);
%title('Absolute error for Chebyschev vs. spline interpolation')

