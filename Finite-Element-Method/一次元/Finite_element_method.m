function [condA,u,y]=Finite_element_method(n,x)
% n 为基函数个数，N 为样点个数
%%
%初始化
m=100000;   N=length(x);
hh=1/n;

h=zeros(n,1);
A=zeros(n,n);
b=zeros(n,1);
p=zeros(n+1,1);
e=zeros(m+1,1);
w=zeros(m+1,1);
x1=zeros(n,1);
x2=zeros(n,1);
u=zeros(N,1);


%%
%网络剖分
for i=1:1:(n+1)
    p(i)=(i-1)*hh;
end


%%
%构造基函数
%我直接构造在u(x)中了


%%
%形成有限元方程
H=1.0/m;
for k=2:2:m
    w(k)=4*H/3;
    w(k+1)=2*H/3;
end
w(1)=H/3; w(m+1)=H/3;

for i=1:1:(n-1)
    x1(i)=p(i);
    x2(i)=p(i+1);
    h(i)=x2(i)-x1(i);
    a1=0; a2=0; a3=0; a4=0; b1=0; b2=0;

    for k=1:1:m+1
        e(k)=(k-1)/m;
        a1=a1+w(k)*(-n+h(i)*pi*pi*(1-e(k))*e(k)/4);
        a2=a1;
        a3=a3+w(k)*(n+h(i)*pi*pi*e(k)*e(k)/4);
        a4=a4+w(k)*(n+h(i)*pi*pi*(1-e(k))*(1-e(k))/4);
        b1=b1+w(k)*h(i)*sin(pi*(x2(i)+h(i)*e(k))/2)*(1-e(k))*pi*pi/2;
        b2=b2+w(k)*h(i)*sin(pi*(x2(i)+h(i)*e(k))/2)*e(k)*pi*pi/2;

    end
    %求解A b
    A(i,i+1)=A(i,i+1)+a1;
    A(i+1,i)=A(i+1,i)+a2;
    b(i)=b(i)+b1;

    A(i,i)=A(i,i)+a3;
    A(i+1,i+1)=A(i+1,i+1)+a4;
    b(i+1)=b(i+1)+b2;

end

h(n)=1/n;
for k=1:1:m+1
    A(1,1)=A(1,1)+w(k)*(n+h(1)*(1-e(k))*(1-e(k))*pi*pi/4);
    b(1)=b(1)+w(k)*h(1)*sin(pi*(x1(1)+h(1)*e(k))/2)*e(k)*pi*pi/2;
end

%%
%求解系数矩阵
y=A\b;
condA = cond(A,2);


%%
%求解有限元解

for i=1:1:N
    ii=fix(x(i)*n)+1;
    if ii==1
        u(i)=y(ii)*(1+n*(x(i)-p(ii+1)));
    end
    if ii>1
        u(i)=y(ii-1)*(1-n*(x(i)-p(ii)))+y(ii)*(1+n*(x(i)-p(ii+1)));
    end
end

end
