## **Code 8.9.2**

==PB21010452 肖羿==

### **1.目标**

使用 Matlab 实现有限差分法，分别求解线性两点边值问题 $\begin{cases}x''=-x\\x(0)=3,x(\frac{\pi}{2}) =7\end{cases}$、$\begin{cases}x''=2e^t-x\\x(0)=2,x(1) =e+\cos 1\end{cases}$，并计算误差 $Error=\max\limits_{t}\{|x_{real}(t)-x(t)|\}$ 和收敛阶 $Ord=\frac{\ln{\frac{Error_{old}}{Error_{now}}}}{\ln{\frac{N_{now}}{N_{old}}}}$。$n$ 分别取 $10,20,40,80,160$ 。

### **2.算法**

真实解分别为 $x_{real}(t)=7\sin t+3\cos t$、$x_{real}(t)=e^t+\cos t$ 。

令 $x''=f(t,x,x')=u(t)+v(t)x+w(t)x'$，两点边值为 $x(a)=\alpha,x(b)=\beta$。

取等距节点 $h=\frac{b-a}{n+1},t_i=a+ih,\forall 0\leqslant i\leqslant n+1$。

解线性方程组：$$\left(\begin{array}{cccccc}d_1 & c_1 & & & & \\a_1 & d_2 & c_2 & & & \\& a_2 & d_3 & c_3 & & \\& & \ddots & \ddots & \ddots & \\& & & a_{n-2} & d_{n-1} & c_{n-1} \\& & & & a_{n-1} & d_n\end{array}\right)\left(\begin{array}{c}x_1 \\x_2 \\x_3 \\\vdots \\x_{n-1} \\x_n\end{array}\right)=\left(\begin{array}{c}b_1-a_0 \alpha \\b_2 \\b_3 \\\vdots \\b_{n-1} \\b_n-c_n \beta\end{array}\right)$$

其中 $u_i=u(t_i),v=v(t_i),w=w(t_i)$，$b_i=-h^2u_i,$ $d_i=2+h^2v_i,$ $a_i=-1-\frac{h}{2}w_{i+1},$ $c_i=-1+\frac{h}{2}w_i$ 。

<div STYLE="page-break-after: always;"></div>
### **3.结果**

|  n   |   (a).error    | (a).order |   (b).error    | (b).order |
| :--: | :------------: | :-------: | :------------: | :-------: |
| 10 | 0.00573244 | -Inf | 0.00029584 | -Inf |
| 20 | 0.00142882 | 2.00433 | 7.39174e-05 | 2.00083 |
| 40 | 0.000357489 | 1.99885 | 1.84965e-05 | 1.99866 |
| 80 | 8.93556e-05 | 2.00027 | 4.62443e-06 | 1.9999 |
| 160 | 2.23407e-05 | 1.99988 | 1.15622e-06 | 1.99986 |

### **4.代码**

```python
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
```