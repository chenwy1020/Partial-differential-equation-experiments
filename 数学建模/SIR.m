clear
clc

n=10000;
h=50/n;
i=zeros(n+1,1); i(1)=0.1;
s=zeros(n+1,1); s(1)=0.9;
r=zeros(n+1,1); r(1)=0;
a=0.4;
b=0.4;


for k=2:1:n+1
    i(k)=i(k-1)+(a*s(k-1)*i(k-1)-b*i(k-1))*h;
    s(k)=s(k-1)+(-a*s(k-1)*i(k-1))*h;
    r(k)=1-i(k-1)-s(k-1);
end


p=(0:1:n)*h;
xlim([0 50]);
ylim([0 1]);
figure
plot(p,i,'DisplayName','i(t)');
hold on
plot(p,s,'DisplayName','s(t)');
hold on
plot(p,r,'DisplayName','r(t)');
hold off
title('a=0.6,b=0.3,s(0)=0.9,i(0)=0.1');
legend("show")




