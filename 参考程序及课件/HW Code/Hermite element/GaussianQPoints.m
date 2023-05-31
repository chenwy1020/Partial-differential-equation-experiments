function [x] = GaussianQPoints(a, b, n)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% 返回高斯型求积公式求积节点
% a, b 为积分区间的端点
% n 为n点高斯型求积公式

if(n==2)
    t=[-sqrt(1/3), sqrt(1/3)];
elseif(n==3)
    t=[-sqrt(3/5), 0, sqrt(3/5)];
elseif(n==4)
    t=[-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7), sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)];
elseif(n==5)
    t=[0.9061798459, 0.5384693101, 0.0000000000 , -0.5384693101,-0.9061798459 ];
end

x=0.5*(a+b+t*(b-a));

end