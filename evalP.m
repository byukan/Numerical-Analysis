function pt=evalP(P,t_int,t_test)
%evaluates the polynomial with coefficients P associated with interpolation
%points t_int at the test points t_test and returns that value in pt.
%task is to evaluate the polynomial P(x)=(x-x0)*Q(x)+b0 at x0.
%work recursively backwards, starting at Fnn(x-x_n-1)
n=length(P)-1;  %gives n, the degree of P
pt=P(n+1)*(t_test-t_int(n));  %compute/starting at Fnn(x-x_n-1)
for j=n:2  %working backwards, using nested multiplication
    pt=(pt+P(j))*(t_test-t_int(j-1));  %compute bj for P
end
pt=pt+P(1);  %output pt, that last term is F(0,0)
end