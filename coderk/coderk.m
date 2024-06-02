RKF54(@(t,x) exp(x*t)+cos(x-t),1,3,0.01,1e-8);
function RKF54(f,t0,x0,h0,delta)
t=t0;x=x0;h=h0;M=50;X=zeros(M,1);Y=X;
fprintf("delta=%g\n",delta);
fprintf("|  i  |   h   |   t   |   x5   |   x4   |   e   |\n");
a1=16/135;a2=0;a3=6656/12825;a4=28561/56430;a5=-9/50;a6=2/55;
b1=25/216;b2=0;b3=1408/2565;b4=2197/4104;b5=-1/5;b6=0;
for i=1:M
    F1=h*f(t,x);
    F2=h*f(t+h/4,x+F1/4);
    F3=h*f(t+h*3/8,x+F1*3/32+F2*9/32);
    F4=h*f(t+h*12/13,x+F1*1932/2197-F2*7200/2197+F3*7296/2197);
    F5=h*f(t+h,x+F1*439/216-F2*8+F3*3680/513-F4*845/4104);
    F6=h*f(t+h/2,x-F1*8/27+F2*2-F3*3544/2565+F4*1859/4104-F5*11/40);
    x5=x+a1*F1+a2*F2+a3*F3+a4*F4+a5*F5+a6*F6;
    x4=x+b1*F1+b2*F2+b3*F3+b4*F4+b5*F5+b6*F6;
    e=abs(x5-x4);t=t+h;x=x5;
    X(i,1)=t;Y(i,1)=x;
    fprintf("| i=%d | %g | %g | %g | %g | %g |\n",i,h,t,x5,x4,e);
    h=0.9*h*(delta/e)^(1/(1+5));
end
figure
plot(X,Y);
xlabel('t')
ylabel('数值解 x(t)')
end
