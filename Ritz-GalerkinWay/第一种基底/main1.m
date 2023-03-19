%������main
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
    c=A\b';%�������Է�����Ľ�
    x=uniform(n,0,1);
    y=assign(n,x,c,N(k));
    err(k)=sqrt(SimpsonCS(n,x,y));%�������

    NcondA(k)=cond(A,2);
    p(k)=(power(k*pi,2)-1)/(power(pi,2)-1);
end

%%��ͼ
plot(N,NcondA,'-o');
xlabel('N');
ylabel('condA');
title('���������� & ϵ������������')

plot(N,err,'-o');
xlabel('N');
ylabel('���err');
title('���������� & ���')

sqrtNcondA=sqrt(NcondA);
plot(N,sqrtNcondA,'-o')   %�����A ���������Ŀ�����N�����Թ�ϵ
xlabel('N');
ylabel('sqrt(condA)');
title('���������� & �������Ŀ���')
