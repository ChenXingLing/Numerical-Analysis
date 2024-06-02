N=8;error_old=0;x_real=-1;
for O=3:8
    x=AdamsBashforth5(@(t,x) (t-exp(-t))/(x+exp(x)),0,0,1,N);
    error=abs(x-x_real);
    fprintf("| %d | %g | %g | %g |\n",N,x,error,log(error_old/error)/log(2));
    error_old=error;N=N*2;
end
function ansx=AdamsBashforth5(f,t0,x0,xn,N)
h=(xn-x0)/N;X14=RungeKutta4(f,t0,x0,h,100);
x=zeros(N+1,1);x(1,1)=x0;x(2:5,1)=X14(1:4,1);
for i=5:N
    x(i+1,1)=x(i,1)+h/720*(1901*f(t0+(i-1)*h,x(i,1))-2774*f(t0+(i-2)*h,x(i-1,1))+2616*f(t0+(i-3)*h,x(i-2,1))-1274*f(t0+(i-4)*h,x(i-3,1))+251*f(t0+(i-5)*h,x(i-4,1)));
end
ansx=x(N+1,1);
end

function X=RungeKutta4(f,t0,x0,h,M)
t=t0;x=x0;X=zeros(4,1);h=h/M;
for i=1:4*M
    F1=h*f(t,x);
    F2=h*f(t+h/2,x+F1/2);
    F3=h*f(t+h/2,x+F2/2);
    F4=h*f(t+h,x+F3);
    x=x+(F1+2*F2+2*F3+F4)/6;
    %fprintf("i=%d: x(%g)=%g\n",i,t,x);
    t=t+h;
    if mod(i,M)==0
        X(round(i/M),1)=x;
    end
end
end

function X=RungeKutta5(f,t0,x0,h,M)
t=t0;x=x0;X=zeros(4,1);h=h/M;
for i=1:4*M
    F1=h*f(t,x);
    F2=h*f(t+h/2,x+F1/2);
    F3=h*f(t+h/2,x+F1/4+F2/4);
    F4=h*f(t+h,x-F2+2*F3);
    F5=h*f(t+2/3*h,x+7/27*F1+10/27*F2+1/27*F4);
    F6=h*f(t+h/5,x+28/625*F1-1/5*F2+546/625*F3+54/625*F4-378/625*F5);
    x=7*x+1/24*F1+5/48*F4+27/56*F5+125/336*F6;
    %fprintf("i=%d: x(%g)=%g\n",i,t,x);
    t=t+h;
    if mod(i,M)==0
        X(round(i/M),1)=x;
    end
end
end