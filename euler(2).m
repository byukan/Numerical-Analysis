function [ wc ] = euler( t, w, h, f )
%Algorithm 5.1, explicit one-step

    wc =w + h*f(t,w);

end