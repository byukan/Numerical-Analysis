function yt = odesolver(a, b, y0, h, f, argf, option)



N = (b-a)/h;
 
 
 U1 = zeros(1,N+1);
 U2 = zeros(1,N+1);
 V1 = zeros(1,N+1);
 V2 = zeros(1,N+1);


 
 wi = [y0,0,0,1];
 U1(1) = wi(1);
 U2(1) = wi(2);
 V1(1) = wi(3);
 V2(1) = wi(4);
 
% StEP 2
 for I = 2 : N+1 
  
  
 wi = rk4(a,wi(1),wi(2),wi(3),wi(4),h,I,@f.dfdy,@f.dfy,@f.r); 

 
 U1(I) = wi(1);
 U2(I) = wi(2);
 V1(I) = wi(3);
 V2(I) = wi(4);
 end;


yt{1,1} = U1;
yt{2,1} = U2;
yt{3,1} = V1;
yt{4,1} = V2;



 
 
 





 function w = rk4(inti,u1,u2,v1,v2,h,i,P,Q,R)
       
      
 x = inti+(I-2)*h;
 t = x+0.5*h;

 K11 = h*u2;
 K12 = h*(P(x)*u2+Q(x)*u1+R(x));
 K21 = h*(u2+0.5*K12);
 K22 = h*(P(t)*(u2+0.5*K12)+Q(t)*(u1+0.5*K11)+R(t));
 K31 = h*(u2+0.5*K22);
 K32 = h*(P(t)*(u2+0.5*K22)+Q(t)*(u1+0.5*K21)+R(t));
 
 t = x+h;
 K41 = h*(u2+K32);
 K42 = h*(P(t)*(u2+K32)+Q(t)*(u1+K31)+R(t));
 u1 = u1+(K11+2*(K21+K31)+K41)/6;
 u2 = u2+(K12+2*(K22+K32)+K42)/6;
 
 K11 = h*v2;
 K12 = h*(P(x)*v2+Q(x)*v1);
 t = x+0.5*h;
 K21 = h*(v2+0.5*K12);
 K22 = h*(P(t)*(v2+0.5*K12)+Q(t)*(v1+0.5*K11));
 K31 = h*(v2+0.5*K22);
 K32 = h*(P(t)*(v2+0.5*K22)+Q(t)*(v1+0.5*K21));
 t = x+h;
 K41 = h*(v2+K32);
 K42 = h*(P(t)*(v2+K32)+Q(t)*(v1+K31));
 
 
 v1 = v1+(K11+2*(K21+K31)+K41)/6;
 v2 = v2+(K12+2*(K22+K32)+K42)/6;
 
 
 w(1) = u1;
 w(2) = u2;
 w(3) = v1;
 w(4) = v2;
 

end






end