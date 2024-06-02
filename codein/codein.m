global n;n=5;
for O=1:8
    int_p1=0;int_p2=0;
    for k=1:n+1
        int_p1=int_p1+f(1-2*(k-1)/n)*integral(@(x) l1(x,k),-1,1);
        int_p2=int_p2+f(-cos(k/(n+2)*pi))*integral(@(x) l2(x,k),-1,1);
    end
    int_f=integral(@(x) f(x),-1,1);
    fprintf("| N=%d | %g %g | %g | %g %g |\n",n,int_p1,int_p2,int_f,abs(int_p1-int_f),abs(int_p2-int_f));
    n=n+5;
end

function Ans_f = f(x)
Ans_f=1./(1+25.*x.*x);
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