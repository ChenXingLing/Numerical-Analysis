## **Code5**

### **1.目标**

使用 Matlab 编写复化 Simpson 积分公式和复化梯形积分公式，分别计算积分 $\int_{0}^{4}\sin x\ dx,$ $\int_{0}^{2\pi}\sin x\ dx$，节点个数取 $N=1,2,4,...,4096$。

### **2.算法**

$$\begin{cases}\text{Simpson}: &\int_{a}^{b}f(x)\ dx=\frac{h}{3}(f(a)+2\sum\limits_{i=2}^{\frac{n}{2}}f(x_{2i-2})+4\sum\limits_{i=1}^{\frac{n}{2}}f(x_{2i-1})+f(b))\\
\text{Gradient}:&\int_{a}^{b}f(x)\ dx=\frac{h}{2}(f(a)+2\sum\limits_{i=1}^{n-1}f(x_i)+f(b))\end{cases}$$

其中 $h=\frac{b-a}{n}$，$x_i=a+ih$。

利用公式计算收敛阶 $Ord=\frac{\ln{\frac{Error_{old}}{Error_{now}}}}{\ln{\frac{N_{now}}{N_{old}}}}$ 。

### **3.结果**

| （积分1） | $\text{Simpson Error}$ | $\text{Simpson Ord}$ | $\text{Gradient Error}$ |$\text{Gradient Ord}$|
| ---- | ---- | ---- | ---- | ---- |
| N=2 | 0.266615 | $-$ | 0.591851 | $-$ |
| N=4 | 0.0104085 | 4.67892 | 0.140156 | 2.07820 |
| N=8 | 0.000591731 | 4.13668 | 0.0345953 | 2.01839 |
| N=16 | 3.61551e-05 | 4.03267 | 0.00862171 | 2.00453 |
| N=32 | 2.24708e-06 | 4.00808 | 0.00215374 | 2.00113 |
| N=64 | 1.40246e-07 | 4.00201 | 0.00053833 | 2.00028 |
| N=128 | 8.76234e-09 | 4.00050 | 0.000134576 | 2.00007 |
| N=256 | 5.47597e-10 | 4.00013 | 3.36436e-05 | 2.00002 |
| N=512 | 3.42233e-11 | 4.00006 | 8.41087e-06 | 2.00000 |
| N=1024 | 2.13962e-12 | 3.99955 | 2.10272e-06 | 2.00000 |
| N=2048 | 1.33449e-13 | 4.00300 | 5.25679e-07 | 2.00000 |
| N=4096 | 5.77316e-15 | 4.53078 | 1.3142e-07 | 2.00000 |

| （积分1） | $\text{Simpson Error}$ | $\text{Simpson Ord}$ | $\text{Gradient Error}$ |$\text{Gradient Ord}$|
| ---- | ---- | ---- | ---- | ---- |
| N=2 | 3.12001e-16 | $-$ | 5.55112e-17 | $-$ |
| N=4 | 5.55112e-17 | 2.49070 | 5.55112e-17 | 0.00000 |
| N=8 | 1.77013e-16 | -1.67301 | 1.18882e-16 | -1.09868 |
| N=16 | 5.55112e-17 | 1.67301 | 1.19128e-17 | 3.31895 |
| N=32 | 6.27775e-17 | -0.17747 | 4.96945e-16 | -5.38250 |
| N=64 | 1.91756e-16 | -1.61095 | 2.97363e-16 | 0.74086 |
| N=128 | 7.41313e-17 | 1.37112 | 4.21669e-16 | -0.50389 |
| N=256 | 3.94307e-16 | -2.41116 | 5.65237e-16 | -0.42274 |
| N=512 | 3.31265e-16 | 0.25133 | 4.1194e-16 | 0.45642 |
| N=1024 | 7.00808e-16 | -1.08103 | 3.14141e-17 | 3.71295 |
| N=2048 | 3.32466e-16 | 1.07581 | 7.96464e-16 | -4.66412 |
| N=4096 | 4.50285e-16 | -0.43763 | 2.0461e-17 | 5.28266 |

### **4.代码**

```python
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

```