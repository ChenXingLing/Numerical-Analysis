eps=0.000001;

fprintf("eps=%.6f\n x      phi(x)        k\n",eps);

%暴力枚举x=0
%ans=pi*pi/6.0;
%n=1;tmp=1.0;err=ans-tmp;
%while err>eps&&n<=300000000
%    n=n+1;
%    tmp=tmp+1.0/n/n;
%    err=err-1.0/n/n;
%end
%fprintf("0.0  %d  %d\n",tmp,n);

%x=0.0~1.0
n=1000000;
for p=0:10
    as=0.0;
    for i=1:n
        as=as+1.0/i/(i+p/10.0);
    end
    fprintf("%.1f  %d  %d\n",p/10.0,as,n);
end

%暴力枚举x=1
%k=1;sum=1.0;error=1.0/(k+1.0);
%while error>=eps
%    k=k+1; error=error-1.0/k+1.0/(k+1);
%end
%fprintf("1.0  %d  %d\n",sum-error,k);

cnt=0;sum=0.0;error=0.0;k=1000000;
for x=1:300
    sum=sum+1.0/x; error=error+1.0/(k+x);
    %fprintf("x=%d k=%d error=%.10f\n",x,k,sum);
    while (error+1.0/k-1.0/(k+x))/x<=eps
        error=error+1.0/k-1.0/(k+x); k=k-1;
    end
    cnt=mod(cnt+1,10);
    if cnt==0
        fprintf("%.1f  %d  %d\n",x,(sum-error)/x,k);
    end
end