## **Codein**

### **1.目标**

使用 Matlab 构造函数 $f(x)=\frac{1}{1+25x^2},x\in[-1,1]$ 的 Lagrange 差值多项式 $p_L(x)$，利用 $\int_{-1}^{1}p_L(x)dx$ 计算积分 $\int _{-1}^{1}f(x)dx$ 的近似值， 并计算误差 $\left| \int_{-1}^1 p_L(x)dx-\int_{-1}^1 f(x)dx \right|$。插值节点取 $x_i=1-\frac{2i}{n}$ 和 $x'_i=-\cos(\frac{i+1}{n+2}\pi)$，$n$ 分别取 $5,10,15,20,25,30,35,40$。

### **2.算法**

构造 $p(x)=\sum\limits_{i=0}^nf(x_i)l_i(x)$，其中 $l_i(x)=\prod\limits_{j=0\\j\neq i}^{n}\frac{x-x_j}{x_i-x_j}$ 。

$\int _{-1}^{1}p_L(x)dx=\sum\limits_{i=0}^{n}A_if(x_i),$ 其中 $A_i=\int _{-1}^{1}l_i(x)dx$。

分别构造拉格朗日基函数 `l1` 和 `l2`，使用 Matlab 内置函数 `integral` 计算积分。

### **3.结果**

| （第一组） | $\int _{-1}^{1}p_L(x)dx$ | $\int _{-1}^{1}f(x)dx$ | $|\int _{-1}^{1}p_L(x)dx-\int _{-1}^{1}f(x)dx|$ |
| :--: | :--: | :--: | :--: |
| N=5 | 0.461538 | 0.54936 | 0.0878218 |
| N=10 | 0.93466 | 0.54936 | 0.3853 |
| N=15 | 0.831112 | 0.54936 | 0.281751 |
| N=20 | -5.36991 | 0.54936 | 5.91927 |
| N=25 | -5.39986 | 0.54936 | 5.94922 |
| N=30 | 153.798 | 0.54936 | 153.249 |
| N=35 | 173.88 | 0.54936 | 173.331 |
| N=40 | -4912.42 | 0.54936 | 4912.97 |

| （第二组） | $\int _{-1}^{1}p_L(x)dx$ | $\int _{-1}^{1}f(x)dx$ | $|\int _{-1}^{1}p_L(x)dx-\int _{-1}^{1}f(x)dx|$|
| :--: | :--: | :--: | :--: |
| N=5 | 0.48114 | 0.54936 | 0.0682199 |
| N=10 | 0.554086 | 0.54936 | 0.00472538 |
| N=15 | 0.547586 | 0.54936 | 0.00177418 |
| N=20 | 0.550011 | 0.54936 | 0.000650485 |
| N=25 | 0.54936 | 0.54936 | 3.94693e-07 |
| N=30 | 0.549322 | 0.54936 | 3.86652e-05 |
| N=35 | 0.549357 | 0.54936 | 3.23414e-06 |
| N=40 | 0.549365 | 0.54936 | 4.50617e-06 |

### **4.代码**

```python
global n;n=5;
for O=1:8
    int_p1=0;int_p2=0;
    for k=1:n+1
        int_p1=int_p1+f(1-2*(k-1)/n)*integral(@(x) l1(x,k),-1,1);
        int_p2=int_p2+f(-cos(k/(n+2)*pi))*integral(@(x) l2(x,k),-1,1);
    end
    int_f=integral(@(x) f(x),-1,1);
    fprintf("| N=%d | %.10f %.10f | %.10f | %g %g |\n",n,int_p1,int_p2,int_f,abs(int_p1-int_f),abs(int_p2-int_f));
    n=n+5;
end

function Ans_f = f(x)
Ans_f=1./(1+x.*x);
end

%拉格朗日基函数1
function Lk=l1(y,k)
global n;Lk=ones(size(y));
for j=1:n+1
    if j~=k
        Lk=Lk.*(y-(1-2*(j-1)/n))./((1-2*(k-1)/n)-(1-2*(j-1)/n));
    end
end
end
%拉格朗日基函数2
function Lk=l2(y,k)
global n;Lk=ones(size(y));
for j=1:n+1
    if j~=k
        Lk=Lk.*(y+cos(j/(n+2)*pi))./(-cos(k/(n+2)*pi)+cos(j/(n+2)*pi));
    end
end
end
```