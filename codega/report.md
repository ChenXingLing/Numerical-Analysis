## **Code-ga**

### **1.目标**

使用 Matlab 编写复化梯形积分公式和复化3点Gauss积分公式，分别计算积分 $\int_0^1e^{-x^2}\ dx$、$\int_0^4\frac{1}{1+x^2}\ dx$、$\int_0^{2\pi}\frac{1}{2+\cos x}\ dx$，节点个数取 $N=1,2,4,...,128$。

### **2.算法**

- $\text{复化 Gradient}:\ \int_{a}^{b}f(x)\ dx=\frac{h}{2}(f(a)+2\sum\limits_{i=1}^{n-1}f(x_i)+f(b))$，其中 $h=\frac{b-a}{n}$，$x_i=a+ih$。

- 做变量代换 $x=\frac{b-a}{2}t+\frac{a+b}{2}$ 得 $\int_a^b f(x)\ dx=\frac{b-a}{2}\int_{-1}^1 f(\frac{b-a}{2}t+\frac{a+b}{2})\ dt$  
再根据三点高斯积分公式 $\int_{-1}^1g(x)\ dx=\sum\limits_{i=0}^{2}A_ig(x_i)=\frac{5}{9}g(-\sqrt{\frac{3}{5}})+\frac{8}{9}g(0)+\frac{5}{9}g(\sqrt{\frac{3}{5}})$ 得：  
$$\text{任意区间 3点 Guass}:\int_a^b f(x)\ dx=\frac{b-a}{2}\left( \frac{5}{9}f(-\sqrt{\frac{3}{5}}\frac{b-a}{2}+\frac{a+b}{2})+\frac{8}{9}f(\frac{a+b}{2})+\frac{5}{9}f(\sqrt{\frac{3}{5}}\frac{b-a}{2}+\frac{a+b}{2}) \right)$$

- 利用公式计算收敛阶 $Ord=\frac{\ln{\frac{Error_{old}}{Error_{now}}}}{\ln{\frac{N_{now}}{N_{old}}}}$ 。

### **3.结果**

| （积分1） | $\text{Gradient Error}$ | $\text{Gradient Ord}$ | $\text{Gauss Error}$ |$\text{Gauss Ord}$|
| :----: | :----: | :----: | :----: | :----: |
| N=2 | 0.0154539 | - |  3.61106e-08 | - |
| N=4 | 0.00384004 | 2.00878 |  4.02152e-10 | 6.48853 |
| N=8 | 0.000958518 | 2.00224 |  5.74218e-12 | 6.13000 |
| N=16 | 0.000239536 | 2.00056 |  8.77076e-14 | 6.03275 |
| N=32 | 5.98782e-05 | 2.00014 |  1.44329e-15 | 5.92527 |
| N=64 | 1.49692e-05 | 2.00004 |  0 | Inf |
| N=128 | 3.74227e-06 | 2.00001 |  2.22045e-16 | -Inf |

<div STYLE="page-break-after: always;"></div>
| （积分2） | $\text{Gradient Error}$ | $\text{Gradient Ord}$ | $\text{Gauss Error}$ |$\text{Gauss Ord}$|
| :----: | :----: | :----: | :----: | :----: |
| N=2 | 0.133006 | - |  0.00012676 | - |
| N=4 | 0.0035941 | 5.20972 |  0.000125931 | 0.00947 |
| N=8 | 0.000564261 | 2.67120 |  2.45799e-07 | 9.00094 |
| N=16 | 0.000144082 | 1.96947 |  2.07057e-12 | 16.85709 |
| N=32 | 3.6038e-05 | 1.99930 |  4.55191e-14 | 5.50741 |
| N=64 | 9.01059e-06 | 1.99982 |  6.66134e-16 | 6.09452 |
| N=128 | 2.25272e-06 | 1.99996 |  2.22045e-16 | 1.58496 |

| （积分3） | $\text{Gradient Error}$ | $\text{Gradient Ord}$ | $\text{Gauss Error}$ |$\text{Gauss Ord}$|
| :----: | :----: | :----: | :----: | :----: |
| N=2 | 0.561191 | - |  0.00611656 | - |
| N=4 | 0.0375927 | 3.89997 |  0.000738328 | 3.05039 |
| N=8 | 0.000192788 | 7.60729 |  4.32607e-06 | 7.41506 |
| N=16 | 5.12258e-09 | 15.19979 |  1.15023e-10 | 15.19885 |
| N=32 | 0 | Inf |  1.33227e-15 | 16.39768 |
| N=64 | 2.22045e-15 | -Inf |  8.88178e-16 | 0.58496 |
| N=128 | 8.88178e-16 | 1.32193 |  1.33227e-15 | -0.58496 |

### **4.代码**

```python
fprintf("I1(f):\n");calc(@(x) exp(-x.*x),0,1);
fprintf("I2(f):\n");calc(@(x) 1./(1+x.*x),0,4);
fprintf("I3(f):\n");calc(@(x) 1./(2+cos(x)),0,2*pi);

function calc(f,a,b)
N=1;int_f=integral(f,a,b);old_GD=0;old_GS=0;
fprintf("             Gradient              Gauss-3 \n");
fprintf("          Error      Ord       Error      Ord\n");
for k=1:7
    N=N*2;
    Grad=Gradient(f,a,b,N);error_GD=abs(int_f-Grad);
    Gau3=Gauss_3(f,a,b,N);error_GS=abs(int_f-Gau3);
    fprintf("| N=%d | %g | %.5f | ",N,error_GD,log(old_GD/error_GD)/log(2));
    fprintf(" %g | %.5f |\n",error_GS,log(old_GS/error_GS)/log(2));
    old_GD=error_GD;old_GS=error_GS;
end
fprintf("\n");
end

function L=Gauss_3_(f,a,b)
x=(b-a)/2*[-sqrt(0.6) 0 sqrt(0.6)]+(a+b)/2;
L=(b-a)/2*f(x)*[5/9;8/9;5/9];
end

function L=Gauss_3(f,a,b,n)
h=(b-a)/n;L=0;
for i=0:n-1
    L=L+Gauss_3_(f,a+i*h,a+(i+1)*h);
end
end

function L=Gradient(f,a,b,n)
h=(b-a)/n;L=f(a)+f(b);
for i=1:n-1
    L=L+2*f(a+i*h);
end
L=L*h/2;
end

```