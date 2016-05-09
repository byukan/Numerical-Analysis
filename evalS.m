function st=evalS(S,t_int,t_test)
%one method to evaluate the spline with coefficients stored in S at the
%point(s) t_test and returns that value in st.
n=length(t_test)-1;
st=0*t_int;
x=t_test;  xx=t_int;  %declare the test points and interpolation points
a=S(:,1);  %fist column of S gets assigned to a
b=S(:,2);  %second column of S gets assigned to b
c=S(:,3);  %third column of S gets assigned to c
d=S(:,4);  %fourth column of S gets assigned to d
for i=1:n
  ix=xx>=x(i) & xx<=x(i+1);
  st(ix)=a(i)+b(i)*(xx(ix)-x(i))+c(i)*(xx(ix)-x(i)).^2+d(i)*(xx(ix)-x(i)).^3;
end