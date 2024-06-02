n=5;
for O=1:4
    x1=zeros(n+1,1);
    for i=0:n
       x1(i+1,1)=5-10.0*i/n;
    end
    x2=zeros(n+1,1);
    for i=0:n
       x2(i+1,1)=-5*cos((2*i+1)/(2*n+2)*pi);
    end
    error1=zeros(101,1);error2=error1;
    for i=0:100
        y=i/10.0-5;
        error1(i+1,1)=abs(f(y)-calc(y,x1));
        error2(i+1,1)=abs(f(y)-calc(y,x2));
    end
    if n==10
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
Ans_f=1.0/(1+x*x);
end

function py=calc(y,x)
py=0.0;n=size(x,1);
for k=1:n
    Lk=1;
    for j=1:n
        if j~=k
            Lk=Lk*(y-x(j,1))/(x(k,1)-x(j,1));
        end
    end
    py=py+f(x(k,1))*Lk;
end
end