## **Code3**

### **1.目标**

使用 Matlab 对函数 $f(x)=e^{x},x\in[0,1]$ 分别构造等距节点的一次线性样条和满足 $S'(0)=1,S'(1)=e$ 的三次样条，并计算误差 $Error=\max\limits_{i}\{|f(x_{i-\frac{1}{2}})-S(x_{i-\frac{1}{2}})|,i=1,...,n\}$ 和收敛阶 $Ord=\frac{\ln{\frac{Error_{old}}{Error_{now}}}}{\ln{\frac{N_{now}}{N_{old}}}}$。$n$ 分别取 $5,10,20,40$ 。

### **2.算法**

取 $h=\frac{1}{n},x_i=\frac{i}{n},y_i=f(x_i),i=0,1,2,...n$。

构造一次线性样条为 $S_i(x)=a_ix+b_i,i=0,1,...n$，其中 $\begin{cases}a_i=n(y_{i+1}-y_{i})\\b_i=y_{i+1}-a_{i}\frac{i+1}{n}\end{cases}$ 。

根据 $\begin{cases}hM_{i-1}+4hM_i+hM_{i+1}=\frac{6}{h}(y_{i+1}-y_i)-\frac{6}{h}(y_i-y_{i-1}),i=1,2,...n-1\\2M_0+M_1=\frac{6}{h}(\frac{y_1-y_0}{h}-S'(x_0))\\M_{n-1}+2M_n=\frac{6}{h}(S'(x_n)-\frac{y_n-y_{n-1}}{h}) \end{cases}$ 建立方程组，解得 $M_i,i=0,1,...n$ 。

则三次样条为 $S_i(x)=y_i+(x-x_i)(C_i+(x-x_i)(B_i+(x-x_i)A_i))$，其中 $\begin{cases}A_i=\frac{M_{i+1}-M_i}{6h}\\B_i=\frac{M_i}{2}\\C_i=\frac{y_{i+1}-y_i}{h}-\frac{h}{6}M_{i+1}-\frac{h}{3}M_i \end{cases}$ 。

<div STYLE="page-break-after: always;"></div>
### **3.结果**

|  n   | Method (1) error |   order   | Method (2) error |   order   |
| :--: | :--------------: | :-------: | :--------------: | :-------: |
| $5$  |  $0.0123082673$  |    $-$    |  $0.0000109074$  |    $-$    |
| $10$ |  $0.0032328105$  | $1.92877$ |  $0.0000006956$  | $3.97094$ |
| $20$ |  $0.0008285329$  | $1.96416$ |  $0.0000000439$  | $3.98688$ |
| $40$ |  $0.0002097304$  | $1.98202$ |  $0.0000000028$  | $3.99379$ |

### **4.代码**

```python
n=5;
error_1=0;
error_2=0;
for O=1:4
    t=zeros(n+1,1); y=t;
    for i=0:n
        t(i+1,1)=1.0*i/n;
        y(i+1,1)=f(t(i+1,1));
    end
    
    %一次
    a=zeros(n,1); b=zeros(n,1);
    for i=0:n-1
        a(i+1,1)=1.0*(y(i+2,1)-y(i+1,1))*n;
        b(i+1,1)=y(i+2,1)-1.0*a(i+1,1)*(i+1)/n;
    end
    error=0;
    for i=1:n
        x=(i-0.5)/n;
        error=max(error,abs(f(x)-a(i,1)*x-b(i,1)));
    end
    fprintf("| %2d | %.10f |",n,error);
    if O>1
        fprintf(" %.5f | ",log(error_1/error)/log(2));
    else
        fprintf("    -    | ");
    end
    error_1=error;

    %三次
    h=1.0/n; u=4*h; b=t; V=zeros(n+1,1); A=zeros(n+1,n+1);
    for i=0:n-1
       b(i+1,1)=6.0*(y(i+2,1)-y(i+1,1))/h;
    end
    V(1,1)=(b(1,1)-6.0*1)/h; % S'(0)=1
    V(n+1,1)=(6.0*exp(1)-b(n,1))/h; % S'(n)=e
    for i=1:n-1
        V(i+1,1)=b(i+1,1)-b(i,1);
        A(i+1,i+1)=u;
        if i<n-1
            A(i+1,i+2)=h; A(i+2,i+1)=h;
        end
    end
    A(1,1)=2; A(1,2)=1; A(2,1)=h;
    A(n+1,n)=1; A(n+1,n+1)=2; A(n,n+1)=h;
    M=eye(n+1)/A*V;
    a=zeros(n,1); b=zeros(n,1); c=zeros(n,1);
    for i=0:n-1
        a(i+1,1)=(M(i+2,1)-M(i+1,1))/6.0/h;
        b(i+1,1)=M(i+1)/2.0;
        c(i+1,1)=(y(i+2,1)-y(i+1,1))/h-h/6.0*M(i+2,1)-h/3.0*M(i+1,1);
    end
    error=0;
    for i=1:n
        x=(i-0.5)/n;
        error=max(error,abs(y(i,1)+(x-t(i,1))*(c(i,1)+(x-t(i,1))*(b(i,1)+(x-t(i,1))*a(i,1)))-f(x)));
    end
    fprintf(" %.10f |",error);
    if O>1
        fprintf(" %.5f | \n",log(error_2/error)/log(2));
    else
        fprintf("    -    | \n");
    end
    error_2=error;
    n=n*2;
end
function Ans_f = f(x)
Ans_f=exp(x);
end
```