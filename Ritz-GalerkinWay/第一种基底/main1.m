%主函数main
clear all
clc
n=10000;
for k=1:1:20
    N(k)=k;
    for i=1:1:N(k)
        b(i)=bieqution(i);
        for j=1:1:N(k)
            A(i,j)=Aijeqution(i,j);
        end
    end
    c=A\b';%计算线性方程组的解
    x=uniform(n,0,1);
    y=assign(n,x,c,N(k));
    err(k)=sqrt(SimpsonCS(n,x,y));%计算误差

    NcondA(k)=cond(A,2);
    p(k)=(power(k*pi,2)-1)/(power(pi,2)-1);
end

%%绘图
plot(N,NcondA,'-o');
xlabel('N');
ylabel('condA');
title('基函数个数 & 系数矩阵条件数')

plot(N,err,'-o');
xlabel('N');
ylabel('误差err');
title('基函数个数 & 误差')

sqrtNcondA=sqrt(NcondA);
plot(N,sqrtNcondA,'-o')   %特殊的A 的条件数的开方与N成线性关系
xlabel('N');
ylabel('sqrt(condA)');
title('基函数个数 & 条件数的开方')
