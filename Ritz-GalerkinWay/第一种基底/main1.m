%������main
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
    c=A\b';%�������Է�����Ľ�
    x=uniform(n,0,1);
    y=assign(n,x,c,N(k));
    err(k)=sqrt(SimpsonCS(n,x,y));%�������

    NcondA(k)=cond(A,2);
end

%%��ͼ
plot(N,NcondA);
xlabel('N');
ylabel('condA');

plot(N,err);
xlabel('N');
ylabel('���err');

sqrtNcondA=sqrt(NcondA);
plot(N,sqrtNcondA)   %�����A ���������Ŀ�����N�����Թ�ϵ
xlabel('N');
ylabel('sqrt(condA)');
