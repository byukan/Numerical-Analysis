function [ wc ] = rk4( ti, wi, h, f )
%Algorithm 5.2, explicit one-step

    k1 = h.*f(ti,wi);
    k2 = h.*f(ti+h/2,wi+k1./2);
    k3 = h.*f(ti+h/2,wi+k2./2);
    k4 = h.*f(ti+h,wi+k3);
    wc = wi + (k1+2.*k2+2.*k3+k4)./6;


end