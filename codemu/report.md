## **Code mu**

PB21010452 肖羿

### **1.目标**

使用 Matlab 实现 $\text{Adams-Bashforth}$ 公式 ，求解初值问题 $\begin{cases}x'=\frac{t-e^{-t}}{x+e^{x}}\\x(0)=0 \end{cases}$ ，该方程真解由等式 $x^2-t^2+2e^x-2e^{-t}=0$ 给出。

### **2.算法**

当 $t=1$ 时，数值求解等式 $x^2-1+2e^x-\frac{2}{e}=0$ 得：$x=-1$ 或 $x\approx -0.155782$ 。

取等距步长 $\bar{h}=\frac{t_N-t_0}{N}=\frac{1}{N}$，$t_i=t_0+i\bar{h}$ $(\forall 1\leqslant i\leqslant N)$ ，其中 $N=8,16,32,64,128,256$ 。

令 $f(t,x)=\frac{t-e^{-t}}{x+e^{x}}$。

五阶 $\text{Adams-Bashforth}$：$x_{n+1}=x_n+\frac{\bar{h}}{720}\left(1901 f_n-2774 f_{n-1}+2616 f_{n-2}\right. 
\left.-1274 f_{n-3}+251 f_{n-4}\right)$，其中 $f_i=f(t_i,x_i)$。

计算前四个初值 $x_1,x_2,x_3,x_4$ 使用四阶 $\text{Runge-Kutta}$：$x(t+h)=x(t)+\frac{1}{6}(F_1+2F_2+2F_3+F_4)$，其中 $\begin{cases}F_1=hf(t,x)\\F_2=hf(t+\frac{h}{2},x+\frac{F_1}{2})\\F_3=hf(t+\frac{h}{2},x+\frac{F2}{2})\\F_4=hf(t+h,x+F_3) \end{cases}$，取步长 $h=\frac{\bar{h}}{100}$ 。

<div STYLE="page-break-after: always;"></div>
### **3.结果**

| n | 数值解 | error | order |
| :--: | :--: | :--: | :--: |
| 8 | -1 | 5.72431e-13 | -Inf |
| 16 | -1 | 1.57216e-11 | -4.77951 |
| 32 | -1 | 2.58593e-12 | 2.604 |
| 64 | -1 | 1.9984e-15 | 10.3376 |
| 128 | -1 | 2.22045e-16 | 3.16993 |
| 256 | -1 | 0 | Inf |

<div STYLE="page-break-after: always;"></div>
### **4.代码**

```python
N=8;error_old=0;x_real=-1;
for O=3:8
    x=AdamsBashforth5(@(t,x) (t-exp(-t))/(x+exp(x)),0,0,1,N);
    error=abs(x-x_real);
    fprintf("| %d | %g | %g | %g |\n",N,x,error,log(error_old/error)/log(2));
    error_old=error;N=N*2;
end
function ansx=AdamsBashforth5(f,t0,x0,xn,N)
h=(xn-x0)/N;X14=RungeKutta4(f,t0,x0,h,100);
x=zeros(N+1,1);x(1,1)=x0;x(2:5,1)=X14(1:4,1);
for i=5:N
    x(i+1,1)=x(i,1)+h/720*(1901*f(t0+(i-1)*h,x(i,1))-2774*f(t0+(i-2)*h,x(i-1,1))+2616*f(t0+(i-3)*h,x(i-2,1))-1274*f(t0+(i-4)*h,x(i-3,1))+251*f(t0+(i-5)*h,x(i-4,1)));
end
ansx=x(N+1,1);
end

function X=RungeKutta4(f,t0,x0,h,M)
t=t0;x=x0;X=zeros(4,1);h=h/M;
for i=1:4*M
    F1=h*f(t,x);
    F2=h*f(t+h/2,x+F1/2);
    F3=h*f(t+h/2,x+F2/2);
    F4=h*f(t+h,x+F3);
    x=x+(F1+2*F2+2*F3+F4)/6;
    %fprintf("i=%d: x(%g)=%g\n",i,t,x);
    t=t+h;
    if mod(i,M)==0
        X(round(i/M),1)=x;
    end
end
end
```