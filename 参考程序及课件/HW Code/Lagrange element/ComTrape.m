function [I] = ComTrape(a, b, f)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

%% 复化梯形公式
% 积分区间 [a, b] 被等分为n份, xj=a+jh, h=(b-a)/n
% length(f)=n+1, f(j)=f(xj)

N=length(f);
h=(b-a)/(N-1);

I=sum(f, 2);
I=I-0.5*f(1)-0.5*f(N);

I=h*I;
end