function [I] = ComSimpsonX(f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% ����Simpson��ʽ
% �������� [a, b] ���ȷ�Ϊn=2m��, xj=a+jh, h=(b-a)/n
% length(f)=n+1, f(j)=f(xj)

N=length(f)-1;
h=1/N;

y=2*ones(N+1,1);
y(2:2:N,1)=4;
y(1, 1)=1;
y(N+1,1)=1;

I=y'*f*h/3;

end