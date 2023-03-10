%主函数main
clear all
clc
n=1000;
for k=1:1:100
    N(k)=2*k;
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
end

%%绘图
plot(N,NcondA);
xlabel('N');
ylabel('condA');

plot(N,err);
xlabel('N');
ylabel('误差err');

sqrtNcondA=sqrt(NcondA);
plot(N,sqrtNcondA)   %特殊的A 的条件数的开方与N成线性关系
xlabel('N');
ylabel('sqrt(condA)');
