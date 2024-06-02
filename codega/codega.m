fprintf("I1(f):\n");calc(@(x) exp(-x.*x),0,1);
fprintf("I2(f):\n");calc(@(x) 1./(1+x.*x),0,4);
fprintf("I3(f):\n");calc(@(x) 1./(2+cos(x)),0,2*pi);

function calc(f,a,b)
N=1;int_f=integral(f,a,b);old_GD=0;old_GS=0;
fprintf("             Gradient              Gauss-3 \n");
fprintf("          Error      Ord       Error      Ord\n");
for k=1:7
    N=N*2;
    Grad=Gradient(f,a,b,N);error_GD=abs(int_f-Grad);
    Gau3=Gauss_3(f,a,b,N);error_GS=abs(int_f-Gau3);
    fprintf("| N=%d | %g | %.5f | ",N,error_GD,log(old_GD/error_GD)/log(2));
    fprintf(" %g | %.5f |\n",error_GS,log(old_GS/error_GS)/log(2));
    old_GD=error_GD;old_GS=error_GS;
end
fprintf("\n");
end

function L=Gauss_3_(f,a,b)
x=(b-a)/2*[-sqrt(0.6) 0 sqrt(0.6)]+(a+b)/2;
L=(b-a)/2*f(x)*[5/9;8/9;5/9];
end

function L=Gauss_3(f,a,b,n)
h=(b-a)/n;L=0;
for i=0:n-1
    L=L+Gauss_3_(f,a+i*h,a+(i+1)*h);
end
end

function L=Gradient(f,a,b,n)
h=(b-a)/n;L=f(a)+f(b);
for i=1:n-1
    L=L+2*f(a+i*h);
end
L=L*h/2;
end
