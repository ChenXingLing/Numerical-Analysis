## **Code2**

### **1.目标**

使用 Matlab 构造函数 $f(x)=\frac{1}{1+25x^2},x\in[-1,1]$ 的 牛顿插值多项式，并计算误差 $\max\limits_{i}\{|f(y_i)-p(y_i)|,y_i=\frac{i}{50}-1,i=0,1,...,100\}$。插值节点取 $x_i=1-\frac{2i}{n}$ 和 $x'_i=-\cos(\frac{2i+1}{2n+2}\pi)$，$n$ 分别取 $5,10,20,40$。

### **2.算法**

构造 $p(x)=\sum\limits_{j=0}^nc_{0j}\prod\limits_{i=0}^{j-1}{x-x_i}$，其中 $c_{ij}=\begin{cases}f(x_i)&,j=0\\\frac{c_{i+1,j-1}-c_{i,j-1}}{x_{i+j}-x_{i}} &,j\neq 0,i+j\leqslant n\end{cases}$ $(0\leqslant i \leqslant n)$

令 $d_i=f(x_i),1\leqslant i \leqslant n$ 为差商表的第一列，现计算第二列并放在 $d_1,d_2,...,d_n$ 的位置，以此类推，最终可得 $d_i=c_{0i},1\leqslant  i\leqslant n$。

### **3.结果**

```python
N=5
Max Error of grid (1) :0.432692
Max Error of grid (2) :0.555911
N=10
Max Error of grid (1) :1.915643
Max Error of grid (2) :0.108929
N=20
Max Error of grid (1) :58.278125
Max Error of grid (2) :0.015325
N=40
Max Error of grid (1) :78689.037829
Max Error of grid (2) :0.000274
```

![](_1.png)

### **4.代码**

```python
n=5;
for O=1:4
    x1=zeros(n+1,1);
    for i=0:n
       x1(i+1,1)=1-2.0*i/n;
    end
    x2=zeros(n+1,1);
    for i=0:n
       x2(i+1,1)=-cos((2*i+1)/(2*n+2)*pi);
    end
    error1=zeros(101,1);error2=error1;
    for i=0:100
        y=i/50.0-1;
        error1(i+1,1)=abs(f(y)-calc(y,x1));
        error2(i+1,1)=abs(f(y)-calc(y,x2));
    end
    if n==20
        x=linspace(-5,5);fy=x;p1y=x;p2y=x;
        m=size(x,2);
        for i=1:m
            fy(1,i)=f(x(1,i));
            p1y(1,i)=calc(x(1,i),x1);
            p2y(1,i)=calc(x(1,i),x2);
        end
        figure
        plot(x,fy,x,p1y,x,p2y);
    end
    fprintf("N=%d\nMax Error of grid (1) :%.6f\nMax Error of grid (2) :%.6f\n",n,max(error1),max(error2));
    n=n*2;
end

function Ans_f = f(x)
Ans_f=1.0/(1+25*x*x);
end

function pX=calc(X,x)
n=size(x,1);
d=zeros(n,1);
for i=1:n
    d(i,1)=f(x(i,1));
end
for j=2:n %计算差商表
    for i=n:-1:j
        d(i,1)=(d(i,1)-d(i-1,1))/(x(i,1)-x(i-j+1,1));
    end
end
pX=d(n,1); %计算p(X)
for i=n-1:-1:1
    pX=d(i,1)+(X-x(i,1))*pX;
end
end
```