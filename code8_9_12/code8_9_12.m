fprintf("          【8.9.2.(a)】       【8.9.2.(b)】\n")
fprintf("   N     error      ord      error      ord\n")
N=10;error_1=0;error_2=0;
for O=1:5
    XY=FDM(@(t) 0,@(t) -1,@(t) 0, 0,3, pi/2,7,N);
    x_real=@(t) 7*sin(t)+3*cos(t);
    %output(N,XY(:,1),x_real(XY(:,1)),XY(:,2));
    error=max(abs(x_real(XY(:,1))-XY(:,2)));
    fprintf("| %g | %g | %g ",N,error,log(error_1/error)/log(2));
    error_1=error;
    XY=FDM(@(t) 2*exp(t),@(t) -1,@(t) 0, 0,2, 1,exp(1)+cos(1),N);
    x_real=@(t) exp(t)+cos(t);
    error=max(abs(x_real(XY(:,1))-XY(:,2)));
    fprintf("| %g | %g |\n",error,log(error_2/error)/log(2));
    error_2=error; 
    N=N*2;
end

function output(N,t,xr,xt)
    fprintf("N=%g:\n",N);
    fprintf("t x_real(t) x(t) error\n");
    for i=1:N+1
        fprintf("| %g | %g | %g | %g |\n",t(i,1),xr(i,1),xt(i,1),abs(xr(i,1)-xt(i,1)));
    end
end
function XY=FDM(u,v,w,x0,y0,xn,yn,N) %有限差分法
XY=zeros(N+1,2);A=zeros(N-1,N-1);B=zeros(N-1,1);
h=(xn-x0)/N;
XY(1,1)=x0;XY(N+1,1)=xn;
XY(1,2)=y0;XY(N+1,2)=yn;
for i=1:N-1
    XY(i+1,1)=x0+h*i;
end
for i=1:N-1
    B(i,1)=-h*h*u(XY(i+1,1)); %b(i)
    A(i,i)=2+h*h*v(XY(i+1,1)); %d(i)
    if i<N-1
        A(i+1,i)=-1-0.5*h*w(XY(i+2,1));%a(i)
        A(i,i+1)=-1+0.5*h*w(XY(i+1,1));%c(i)
    end
end
B(1,1)=B(1,1)-y0*(-1-0.5*h*w(XY(1,1))); %b(1)-a(0)*y0
B(N-1,1)=B(N-1,1)-yn*(-1+0.5*h*w(XY(N,1)));%b(n-1)-c(n-1)*yn
Y=A\B;
for i=1:N-1
    XY(i+1,2)=Y(i,1);
end
end