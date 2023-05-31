function [I] = QuadratureRules(a, b, f, A)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% 区间[a, b] 上的数值积分公式
%% f 为其在求积节点上的取值
%% A 为积分公式类型

if (A==0) % 四点 Gauss 型求积公式
    I=GaussianQuadrature(a, b, f, 4);
elseif (A==1) % 复化 Simpson 公式
    I=ComSimpson(a, b, f);
elseif (A==2) % 复化梯形公式
    I=ComTrape(a, b, f);
end