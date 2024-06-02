fprintf("I1(f):\n");calc(@(x) sin(x)./x,0,1);
fprintf("I2(f):\n");calc(@(x) (cos(x)-exp(x))./sin(x),-1,1);
fprintf("I3(f):\n");calc(@(x) 1./(x.*exp(1./x)),0,1);

function calc(f,a,b)
M=7;int_f=integral(f,a,b);
fprintf("精确值: %.6f\n",int_f);
R=Romberg(f,a,b,M,1e-20);
fprintf("Romberg阵列:\n");
for n=0:M
    for m=0:n
        fprintf("| %.6f  ",R(n+1,m+1));
    end
    fprintf("|\n");
end
fprintf("误差：\n");
for m=1:M
    fprintf("N=%d: error=%g\n",m,abs(R(m+1,m+1)-int_f));
end
fprintf("\n");
end

function R=Romberg(f,a,b,M,eps)
h=b-a;R=zeros(M+1,M+1);
%syms px;
%lim_a=limit(f(px),px,a);
R(1,1)=h/2*[f(a+eps)+f(b-eps)];%每个点向附近偏移eps后进行计算
mi=1;
for n=1:M
    h=h/2;
    R(n+1,1)=R(n,1)/2;
    for i=1:mi
        %calc_f=f(a+(2*i-1)*h);
        %calc_f=limit(f(px),px,a+(2*i-1)*h);
        calc_f=f(a+(2*i-1)*h+eps);
        R(n+1,1)=R(n+1,1)+h*calc_f;
    end
    p=1;
    for m=1:n
        p=p*4;
        R(n+1,m+1)=R(n+1,m)+(R(n+1,m)-R(n,m))/(p-1); 
    end
    mi=mi*2;
end
end
