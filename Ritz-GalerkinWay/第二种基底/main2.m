%主函数main
clear all
clc
n=10000;
for k=1:1:10
    N(k)=k;
    for i=1:1:N(k)
        b(i)=bieqution(i);
        for j=1:1:N(k)
            A(i,j)=Aijeqution(i,j);
        end
    end
    c=A\b';%计算线性方程组的解
    NcondA(k)=cond(A,2);

    x=uniform(n,0,1);
    y=assign(n,x,c,N(k));
    err(k)=sqrt(SimpsonCS(n,x,y));%计算误差

    
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

sqrtNcondA = log(NcondA);
plot(N,sqrtNcondA,'-o')   %特殊的A 的条件数的对数与N基本成线性关系
xlabel('N');
ylabel('log(condA)');
title('基函数个数 & 条件数的对数 ')