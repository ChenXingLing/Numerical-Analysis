calc(0,4);calc(0,2*pi)
function calc(a,b)
N=1;int_f=integral(@sin,a,b);old_S=0;old_G=0;
fprintf("           Simpson               Gradient \n");
fprintf("       Error      Ord        Error       Ord\n");
for k=1:12
    N=N*2;
    Simp=Simpson(@sin,a,b,N);error_S=abs(int_f-Simp);    
    Grad=Gradient(@sin,a,b,N);error_G=abs(int_f-Grad);
    fprintf("N=%d | %g | %.5f |",N,error_S,log(old_S/error_S)/log(2));
    fprintf(" %g | %.5f |\n",error_G,log(old_G/error_G)/log(2));
    old_S=error_S;old_G=error_G;
end
fprintf("\n");
end

function L=Simpson(f,a,b,n)
h=(b-a)/n;x=@(i) a+i*h;
L=f(x(0))+f(x(n));
for i=2:n/2
    L=L+2*f(x(2*i-2));
end
for i=1:n/2
    L=L+4*f(x(2*i-1));
end
L=L*h/3;
end

function L=Gradient(f,a,b,n)
h=(b-a)/n;L=f(a)+f(b);
for i=1:n-1
    L=L+2*f(a+i*h);
end
L=L*h/2;
end
