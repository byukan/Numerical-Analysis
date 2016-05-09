function [ yt ] = hw2odesolver( a, b, y0, hmin, hmax, tol, f, argf, option)
%routine implements an adaptive step-size method ex. rkf45
%Runge-Kutta-Fehlberg algorithm approximates soln of IVP with local
%truncation error within a given tolerance
%part of Algorithm 5.3 in Burden & Faires 

%yt is a Nx3 matrix

switch option
    case 'rkf';
        
    %argPhi is the 8*7 Butcher Tableau matrix for rkf45   
argPhi=[0 0 0 0 0 0 0
    1/4 1/4 0 0 0 0 0
    3/8 3/32 9/32 0 0 0 0
    12/13 1932/2197 -7200/2197 7296/2197 0 0 0
    1 439/216 -8 3680/513 -845/4104 0 0
    1/2 -8/27 2 -3544/2565 1859/4104 -11/40 0
    0 25/216 0 1408/2565 2197/4104 -1/5 0
    0 16/135 0 6656/12825 28561/56430 -9/50 2/55]; 

ti(1)=a;  %set t=a
wi(1,:)=y0;  %y0 is initial condition alpha, vectorized
h=hmax;
j=2;
i=1;  %set up a flag index starting at 1

while i==1  %while i=1, implement the rkf method
    %use (t,w) in next while loop step
    
[wc] = rkf(ti(j-1),wi(j-1),h,f,argf,argPhi);  %wc is a vector where wc(1) is w4
%and wc(2) is w5
  
  %R is the difference between w tilde_i+1 and w_i+1 over h
  R=1/h*abs(wc(:,1)-wc(:,2));
  
  % get L norm (norm infinity) of R, since R is a matrix and we need to
  % compare it to a scalar tolerance
  R = norm(R, inf)
  
  if R<=tol
      ti(j)=ti(j-1)+h;  %approximation accepted
      wi(j,:)=wc(:,1); %w approximates y(t)
      j=j+1;
  end
  delta=0.84*(tol/R)^(1/4);
  if delta<=0.1
    h=0.1*h;
  elseif delta>=4
    h=4*h;
  else
    h=delta*h;  %calculate new h
  end
  if h>hmax, h=hmax; end
  if ti(j-1)>=b
    break;  %exit while loop set i=0, flag=0
  elseif ti(j-1)+h>b
    h=b-ti(j-1);
  elseif h<hmin  %procedure completed unsuccessfully
    error('Minimum h exceeded');
  end
end

yt = [ti' wi'];

end

