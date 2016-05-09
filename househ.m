function a = househ(a)
% Turn matrix a into upper- hessenberg with same eigenvalues
[m,n]=size(a);
if m~=n, error('Not a square matrix!'); end
for i=1:n-2
        c=zeros(n,1);
        c(i+1:n)=a(i+1:m,i); % set up c
        q=norm(c);
        b=zeros(n,1);
        if a(i+1,i)>=0,
            b(i+1)=-q; 
        else
            b(i+1)=q;
        end
        w=(c-b)/norm(c-b);
        Q=eye(n)-2*w*w';
        a=Q*a*Q;
        a(i+2:n,i)=0; %let's correct the roundoff error
end