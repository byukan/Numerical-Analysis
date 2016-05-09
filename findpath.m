function [path, coeff]=findpath(ip,tp,method,varargin)
%finds the path matrix of the form (1) and the coefficient matrix of the
%form (2) or (4).  Here ip is a matrix of the form (1) storing all
%interpolation points, tp is a column vector storing t0,...tm-1, and method
%specifies the interpolation method, eg. strings "lagrange" and "spline"

if strcmp(method,'lagrange')== 1;  %specify lagrange for interpolation method
    %construct & evaluate a single polynomial that interpolates all n nodes
    %first find coefficients of lagrange polynomial that gives x coordinates
    
     [Px,Qx]=divdif(ip(:,1),ip(:,2));  %use columns 1(time) and 2 x(t) of ip, which stores interpolation points
         x_coord=evalP(Px,ip(:,1),tp);  %Px gives x coordinate for each time t

     [Py,Qy]=divdif(ip(:,1),ip(:,3));  %use columns 1(time) and 3 y(t) of ip, which stores interpolation points
         y_coord=evalP(Py,ip(:,1),tp);  %Py gives y coordinate for each time t

            path=[ip(:,1),x_coord,y_coord];  %path matrix that gives the x and y coordinates at each time t
            coeff=[ip(:,1),Px,Py];  %save the polynomial coefficients
end





if strcmp(method,'spline')== 1;  %specify spline for interpolation method
   Sx=spline3(ip(:,1),ip(:,2),varargin); 
   %use spline3 to evaluate Sx using given interpolation points t, and x(t)
   Sy=spline3(ip(:,1),ip(:,3),varargin);  
   %use spline3 to evaluate Sy using given interpolation points t, and y(t)
   
   sx_coord=evalS(Sx,ip(:,1),tp);  %evaluate spline w/ coefficients stored in S
   sy_coord=evalS(Sy,ip(:,1),tp);  %at given points tp, and returns sx,sy coordinates
   
   path=[ip(:,1), sx_coord, sy_coord]; %path matrix that gives x,y coords at each time t
   coeff=[ip(:,1),Sx(:,1),Sx(:,2),Sx(:,3),Sx(:,4),Sy(:,1),Sy(:,2),Sy(:,3),Sy(:,4) ];
   %final matrix of the form (4), save all coefficients
   
end