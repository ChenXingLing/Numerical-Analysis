n=5;
error_1=0;
error_2=0;
for O=1:4
    t=zeros(n+1,1); y=t;
    for i=0:n
        t(i+1,1)=1.0*i/n;
        y(i+1,1)=f(t(i+1,1));
    end
    
    %一次
    a=zeros(n,1); b=zeros(n,1);
    for i=0:n-1
        a(i+1,1)=1.0*(y(i+2,1)-y(i+1,1))*n;
        b(i+1,1)=y(i+2,1)-1.0*a(i+1,1)*(i+1)/n;
    end
    error=0;
    for i=1:n
        x=(i-0.5)/n;
        error=max(error,abs(f(x)-a(i,1)*x-b(i,1)));
    end
    fprintf("| %2d | %.10f |",n,error);
    if O>1
        fprintf(" %.5f | ",log(error_1/error)/log(2));
    else
        fprintf("    -    | ");
    end
    error_1=error;

    %三次
    h=1.0/n; u=4*h; b=t; V=zeros(n+1,1); A=zeros(n+1,n+1);
    for i=0:n-1
       b(i+1,1)=6.0*(y(i+2,1)-y(i+1,1))/h;
    end
    V(1,1)=(b(1,1)-6.0*1)/h; % S'(0)=1
    V(n+1,1)=(6.0*exp(1)-b(n,1))/h; % S'(n)=e
    for i=1:n-1
        V(i+1,1)=b(i+1,1)-b(i,1);
        A(i+1,i+1)=u;
        if i<n-1
            A(i+1,i+2)=h; A(i+2,i+1)=h;
        end
    end
    A(1,1)=2; A(1,2)=1; A(2,1)=h;
    A(n+1,n)=1; A(n+1,n+1)=2; A(n,n+1)=h;
    M=eye(n+1)/A*V;
    a=zeros(n,1); b=zeros(n,1); c=zeros(n,1);
    for i=0:n-1
        a(i+1,1)=(M(i+2,1)-M(i+1,1))/6.0/h;
        b(i+1,1)=M(i+1)/2.0;
        c(i+1,1)=(y(i+2,1)-y(i+1,1))/h-h/6.0*M(i+2,1)-h/3.0*M(i+1,1);
    end
    error=0;
    for i=1:n
        x=(i-0.5)/n;
        error=max(error,abs(y(i,1)+(x-t(i,1))*(c(i,1)+(x-t(i,1))*(b(i,1)+(x-t(i,1))*a(i,1)))-f(x)));
    end
    fprintf(" %.10f |",error);
    if O>1
        fprintf(" %.5f | \n",log(error_2/error)/log(2));
    else
        fprintf("    -    | \n");
    end
    error_2=error;
    n=n*2;
end
function Ans_f = f(x)
Ans_f=exp(x);
end


function Ans_g = g(x,n,a,b)
    for i=0:n-1
        if x>=1.0*i/n&&x<=1.0*(i+1)/n
            Ans_g=a(i+1,1)*x+b(i+1,1);
        end
    end
end