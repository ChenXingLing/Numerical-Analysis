format long
Richardson(@log,1,3,3);
Richardson(@tan,1,asin(0.8),4);
Richardson(@f3,1,0,5);

function Ans_f = Richardson(f,h,x,M)
    D=zeros(M+1,M+1);    
    for i=0:M
        D(i+1,0+1)=(f(x+h)-f(x-h))*0.5/h;
        h=h/2.0;
    end
    for k=1:M
        for n=k:M
            D(n+1,k+1)=((4^k)*D(n+1,k-1+1)-D(n-1+1,k-1+1))/((4^k)-1);
        end
    end
    for i=0:M
        for j=0:i
            fprintf("%.6f ",D(i+1,j+1));
        end
        fprintf("\n");
    end
    Ans_f=D;
end

function Ans_f3 = f3(x)
    Ans_f3=sin(x*x+x/3.0);
end