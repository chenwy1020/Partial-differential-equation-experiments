function [I] = ComSimpson(a, b, f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% 复化Simpson公式
% 积分区间 [a, b] 被等分为n=2m份, xj=a+jh, h=(b-a)/n
% length(f)=n+1, f(j)=f(xj)

N=length(f)-1;
h=(b-a)/N;

y=2*ones(1, N+1);
y(1, 2:2:N)=4;
y(1, 1)=1;
y(1, N+1)=1;

I=f*y'*h/3;

end