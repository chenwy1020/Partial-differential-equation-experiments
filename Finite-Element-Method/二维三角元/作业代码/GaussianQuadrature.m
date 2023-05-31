function [I] = GaussianQuadrature(a, b, f, n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% ��˹�������ʽ
% a, b Ϊ��������Ķ˵�
% f Ϊ��������������ڵ��ϵ�ȡֵ
% n Ϊn���˹�������ʽ

if (n==2)
    I=f(1)+f(2);
elseif (n==3)
    I=5/9*f(1)+8/9*f(2)+5/9*f(3);
elseif(n==4)
   I=(18-sqrt(30))/36*f(1)+(18+sqrt(30))/36*f(2)+(18+sqrt(30))/36*f(3)+(18-sqrt(30))/36*f(4);
elseif(n==5)
    I=0.2369268850*(f(1)+f(5))+0.4786286705*(f(2)+f(4))+0.5688888889*f(3);
end

% ���� Jacobi ����ʽ
I=I*(b-a)/2;
end