n=5;
for O=1:4
    x1=zeros(n+1,1);
    for i=0:n
       x1(i+1,1)=1-2.0*i/n;
    end
    x2=zeros(n+1,1);
    for i=0:n
       x2(i+1,1)=-cos((2*i+1)/(2*n+2)*pi);
    end
    error1=zeros(101,1);error2=error1;
    for i=0:100
        y=i/50.0-1;
        error1(i+1,1)=abs(f(y)-calc(y,x1));
        error2(i+1,1)=abs(f(y)-calc(y,x2));
    end
    if n==20
        x=linspace(-5,5);fy=x;p1y=x;p2y=x;
        m=size(x,2);
        for i=1:m
            fy(1,i)=f(x(1,i));
            p1y(1,i)=calc(x(1,i),x1);
            p2y(1,i)=calc(x(1,i),x2);
        end
        figure
        plot(x,fy,x,p1y,x,p2y);
    end
    fprintf("N=%d\nMax Error of grid (1) :%.6f\nMax Error of grid (2) :%.6f\n",n,max(error1),max(error2));
    n=n*2;
end

function Ans_f = f(x)
Ans_f=1.0/(1+25*x*x);
end

function pX=calc(X,x)
n=size(x,1);
d=zeros(n,1);
for i=1:n
    d(i,1)=f(x(i,1));
end
for j=2:n %计算差商表
    for i=n:-1:j
        d(i,1)=(d(i,1)-d(i-1,1))/(x(i,1)-x(i-j+1,1));
    end
end
pX=d(n,1); %计算p(X)
for i=n-1:-1:1
    pX=d(i,1)+(X-x(i,1))*pX;
end
end