## **Code0**

### **1.目标**

使用 Matlab 计算Hamming级数 $\phi(x)=\sum\limits_{k=1}^{\infty}\frac{1}{k(k+x)}$ 分别在 $x=0.0,0.1,0.2,...,1.0,10.0,20.0,...,300.0$ 时的取值，误差小于 $10^{-6}$，并给出相应的最小 $k$ 值。

### **2.算法**

由 $\phi(0)=\frac{\pi^2}{6},\phi(1)=1$ 可暴力求解相应 $k$ 值均为 $10^6$。

```python
%暴力枚举x=0
ans=pi*pi/6.0;
n=1;tmp=1.0;err=ans-tmp;
while err>eps&&n<=300000000
    n=n+1;
    tmp=tmp+1.0/n/n;
    err=err-1.0/n/n;
end
fprintf("0.0  %d  %d\n",tmp,n);

%暴力枚举x=1
k=1;sum=1.0;error=1.0/(k+1.0);
while error>=eps
    k=k+1; error=error-1.0/k+1.0/(k+1);
end
fprintf("1.0  %d  %d\n",sum-error,k);
```

当 $0<x<1$ 时，误差 $\epsilon(k)= \sum\limits_{i=k+1}^{\infty}\frac{1}{i(i+x)}$ 对于同一 $k$ 值满足 $\epsilon_0(k)>\epsilon_x(k)>\epsilon_1(k)$，所以对应 $k$ 均为 $10^6$。

当 $x>1$ 且为整数时，$\phi(x)=\frac{1}{x}\sum\limits_{i=1}^{x}\frac{1}{i}$，随着 $x$ 的增大，误差 $\epsilon_x(k)=\frac{1}{x}\sum\limits_{i=k+1}^{\infty}\frac{1}{i}$ 减小，最小 $k$ 值减小。

### **3.结果**

![](_1.Png)

![](_2.png)

### **4.代码**

```python
eps=0.000001;

fprintf("eps=%.6f\n x      phi(x)        k\n",eps);

%x=0.0~1.0
n=1000000;
for p=0:10
    as=0.0;
    for i=1:n
        as=as+1.0/i/(i+p/10.0);
    end
    fprintf("%.1f  %d  %d\n",p/10.0,as,n);
end

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
```