# Numerical Analysis

数值分析算法实现（Hamming级数、Lagrange插值、牛顿插值、线性样条、Richardson外推、复化积分、龙贝格、复化3点Gauss、龙格-库塔、RKF、Adams-Bashforth、有限差分）

## **【Code 0】**

使用 Matlab 计算Hamming级数 $\phi(x)=\sum\limits_{k=1}^{\infty}\frac{1}{k(k+x)}$ 分别在 $x=0.0,0.1,0.2,...,1.0,10.0,20.0,...,300.0$ 时的取值，误差小于 $10^{-6}$，并给出相应的最小 $k$ 值。

[【report】](./code0/report.md)

## **【Code 1】**

使用 Matlab 构造函数 $f(x)=\frac{1}{1+x^2},x\in[-5,5]$ 的 Lagrange 差值多项式，并计算误差 $\max\limits_{i}\{|f(y_i)-p(y_i)|,y_i=\frac{i}{10}-5,i=0,1,...,100\}$。插值节点取 $x_i=5-\frac{10i}{n}$ 和 $x'_i=-5\cos(\frac{2i+1}{2n+2}\pi)$，$n$ 分别取 $5,10,20,40$。

[【report】](./code1/report.md)

## **【Code 2】**

使用 Matlab 构造函数 $f(x)=\frac{1}{1+25x^2},x\in[-1,1]$ 的 牛顿插值多项式，并计算误差 $\max\limits_{i}\{|f(y_i)-p(y_i)|,y_i=\frac{i}{50}-1,i=0,1,...,100\}$。插值节点取 $x_i=1-\frac{2i}{n}$ 和 $x'_i=-\cos(\frac{2i+1}{2n+2}\pi)$，$n$ 分别取 $5,10,20,40$。

[【report】](./code2/report.md)

## **【Code 3】**

使用 Matlab 对函数 $f(x)=e^{x},x\in[0,1]$ 分别构造等距节点的一次线性样条和满足 $S'(0)=1,S'(1)=e$ 的三次样条，并计算误差 $Error=\max\limits_{i}\{|f(x_{i-\frac{1}{2}})-S(x_{i-\frac{1}{2}})|,i=1,...,n\}$ 和收敛阶 $Ord=\frac{\ln{\frac{Error_{old}}{Error_{now}}}}{\ln{\frac{N_{now}}{N_{old}}}}$。$n$ 分别取 $5,10,20,40$ 。

[【report】](./code3/report.md)

## **【Code 4】**

使用 Matlab 实现 Richardson 外推算法计算 $f'(x)$（取$h=1$），输出相应三角阵列。

函数 $f(x)$ 分别取：

- $\ln x$（计算 $x=3,M=3$）
- $\tan x$（计算 $x=sin^{-1}(0.8),M=4$）
- $\sin(x^2+\frac{x}{3})$（计算 $x=0,M=5$）

[【report】](./code4/report.md)

## **【Code 5】**

使用 Matlab 编写复化 Simpson 积分公式和复化梯形积分公式，分别计算积分 $\int_{0}^{4}\sin x\ dx,$ $\int_{0}^{2\pi}\sin x\ dx$，节点个数取 $N=1,2,4,...,4096$。

[【report】](./code5/report.md)

## **【Code in】**

使用 Matlab 构造函数 $f(x)=\frac{1}{1+25x^2},x\in[-1,1]$ 的 Lagrange 差值多项式 $p_L(x)$，利用 $\int_{-1}^{1}p_L(x)dx$ 计算积分 $\int _{-1}^{1}f(x)dx$ 的近似值， 并计算误差 $\left| \int_{-1}^1 p_L(x)dx-\int_{-1}^1 f(x)dx \right|$。插值节点取 $x_i=1-\frac{2i}{n}$ 和 $x'_i=-\cos(\frac{i+1}{n+2}\pi)$，$n$ 分别取 $5,10,15,20,25,30,35,40$。

[【report】](./codein/report.md)

## **【Code 7.4】**

使用 Matlab 实现龙贝格算法，分别计算积分 $I_1=\int_0^1\frac{\sin x}{x}\ dx$、$I_2=\int_{-1}^{1}\frac{\cos x-e^x}{\sin x}\ dx$、$I_3=\int_1^{\infty}\frac{1}{xe^x}\ dx$ 并打印龙贝格阵列，计算误差 $error=|I-R(m,m)|,m=1,2,3...,7$ 。

[【report】](./code7_4/report.md)

## **【Code ga】**

使用 Matlab 编写复化梯形积分公式和复化3点Gauss积分公式，分别计算积分 $\int_0^1e^{-x^2}\ dx$、$\int_0^4\frac{1}{1+x^2}\ dx$、$\int_0^{2\pi}\frac{1}{2+\cos x}\ dx$，节点个数取 $N=1,2,4,...,128$。

[【report】](./codega/report.md)

## **【Code 8.3.2】**

使用 Matlab 实现四阶龙格-库塔方法，在区间 $[0,5]$ 上求解初值问题 $\begin{cases}x'=\lambda x+\cos t-\lambda\sin t\\x(0)=0 \end{cases}$，计算数值解和精确解的误差 $error$，步长 $h=0.01$，$\lambda$ 分别取 $5,-5,-10$。

[【report】](./code832/report.md)

## **【Code rk】**

使用 Matlab 实现RKF54方法，求解初值问题 $\begin{cases}x'=e^{xt}+\cos(x-t)\\x(1)=3 \end{cases}$，取初值步长 $h=0.01$。

[【report】](./coderk/report.md)

## **【Code mu】**

使用 Matlab 实现 $\text{Adams-Bashforth}$ 公式 ，求解初值问题 $\begin{cases}x'=\frac{t-e^{-t}}{x+e^{x}}\\x(0)=0 \end{cases}$ ，该方程真解由等式 $x^2-t^2+2e^x-2e^{-t}=0$ 给出。

[【report】](./codemu/report.md)

## **【Code 8.9.2】**

使用 Matlab 实现有限差分法，分别求解线性两点边值问题 $\begin{cases}x''=-x\\x(0)=3,x(\frac{\pi}{2}) =7\end{cases}$、$\begin{cases}x''=2e^t-x\\x(0)=2,x(1) =e+\cos 1\end{cases}$，并计算误差 $Error=\max\limits_{t}\{|x_{real}(t)-x(t)|\}$ 和收敛阶 $Ord=\frac{\ln{\frac{Error_{old}}{Error_{now}}}}{\ln{\frac{N_{now}}{N_{old}}}}$。$n$ 分别取 $10,20,40,80,160$ 。

[【report】](./code8_9_12/report.md)

## **【Code st】**

使用 Mathematica 绘制五阶 Adams-Bashforth 公式和五阶 Adams-Moulton 公式的绝对稳定区域。

[【report】](./codest/report.pdf)