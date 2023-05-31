function [I] = GaussianQuadrature2D1PT(f, n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% 高斯型求积公式
% a, b 为积分区间的端点
% f 为被积函数于求积节点上的取值
% n 为n点高斯型求积公式

if(n==3)
    I=f(1)+f(2)+f(3);
    I=I/3;
elseif(n==4)
   I=25*(f(1)+f(2)+f(3))/48-9*f(4)/16;
elseif(n==6)
    I=f(1)+f(2)+f(3)+f(4)+f(5)+f(6);
    I=I/6;
end

% 乘以 Jacobi 行列式,
% 此时 Jacobi 行列式的值为 S 即三角形面积，放到主程序中进行运算

end