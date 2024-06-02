fprintf("lambda=5:\n");RungeKutta(5);
fprintf("lambda=-5:\n");RungeKutta(-5);
fprintf("lambda=10:\n");RungeKutta(-10);

function RungeKutta(lam)
f=@(t,x) lam*x+cos(t)-lam*sin(t);
real_x=@(t) sin(t);
t=0;x=0;h=0.01;M=500;
for i=1:M
    F1=h*f(t,x);
    F2=h*f(t+h/2,x+F1/2);
    F3=h*f(t+h/2,x+F2/2);
    F4=h*f(t+h,x+F3);
    x=x+(F1+2*F2+2*F3+F4)/6;
    t=t+h;
    error=abs(real_x(t)-x);
    fprintf("| i=%d | %g | %g | %g | %g |\n",i,t,x,real_x(t),error);
end
fprintf("\n");
end
