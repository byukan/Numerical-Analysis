function [wc] = rkf(ti,wi,h,f,argf,argPhi)
%routine that implements the Runge-Kutta-Fehlberg method
%part of Algorithm 5.3 in Burden & Faires

%k 's need to take x and y

  k1=h*f(ti,wi,argf);  %set the K1 to K6 as described in algorithm 5.3
  k2=h*f(ti+argPhi(2,1)*h,wi+k1*argPhi(2,2),argf);  %use argPhi butcher tableau matrix
  k3=h*f(ti+argPhi(3,1)*h,wi+argPhi(3,2)*k1+argPhi(3,3)*k2, argf);
  k4=h*f(ti+argPhi(4,1)*h,wi+argPhi(4,2)*k1+argPhi(4,3)*k2+argPhi(4,4)*k3, argf);
  k5=h*f(ti+h,wi+argPhi(5,2)*k1+argPhi(5,3)*k2+argPhi(5,4)*k3+argPhi(5,5)*k4, argf);
  k6=h*f(ti+h*argPhi(6,1),wi+argPhi(6,2)*k1+argPhi(6,3)*k2+argPhi(6,4)*k3+argPhi(6,5)*k4+argPhi(6,6)*k5, argf);
  
  %Runge-Kutta method with LTE of order five, w5
  wtilde_iplus1 = wi+argPhi(8,2)*k1+argPhi(8,4)*k3+argPhi(8,5)*k4+argPhi(8,6)*k5+argPhi(8,7)*k6;
  %use above to estimate the local error in a Runge-Kutta method of order
  %four givn by:
  w_iplus1 = wi+argPhi(7,2)*k1+argPhi(7,4)*k3+argPhi(7,5)*k4+argPhi(7,6)*k5;
  
  %w4 first row, w5 second row
  wc=[w_iplus1; %the van der pol has 2 variables so we have x and y
      wtilde_iplus1 ]; %and wc is 2x2 matrix
 
  
end