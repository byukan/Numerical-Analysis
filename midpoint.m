        function [] = midpoint
        w1=w(i+1,:)';
        f1=(-1/2)*inv(argf(w1(:)))*F
        w2=w(i+1,:)'+h*f1);
        wc=w(i+1,:)'-h*inv(argf(w2(:)))*F;
        end