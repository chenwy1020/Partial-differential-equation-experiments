function [I] = GaussQuadratureAB(a, b, c, d, f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% 乘积型高斯求积公式：4点公式
% f(s, t)=f(xi_s, \eta_t)

w=[(18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];

I=w*f*w';
I=I*(b-a)/2*(d-c)/2;


end